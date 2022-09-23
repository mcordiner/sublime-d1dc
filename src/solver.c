/*
 *  solver.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
TODO:
  - The test to run calculateJBar() etc in levelPops just tests dens[0]. This is a bit sloppy.
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#include "lime.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

//###
#define Ith(v,i)    NV_Ith_S(v,i)         /* Ith numbers components 0..NEQ-1 */
#define IJth(sunMatrix,i,j) SM_ELEMENT_D(sunMatrix,i,j) /* IJth numbers rows,cols 0..NEQ-1 */

/* Data concerning a single grid vertex which is passed from calculateJBar() to solveStatEq(). This data needs to be thread-safe. */
typedef struct {
  double *jbar,*phot,*vfac,*vfac_loc;
} gridPointData;

struct blend{
  int molJ, lineJ;
  double deltaV;
};

struct lineWithBlends{
  int lineI, numBlends;
  struct blend *blends;
};

struct molWithBlends{
  int molI, numLinesWithBlends;
  struct lineWithBlends *lines;
};

struct blendInfo{
  int numMolsWithBlends;
  struct molWithBlends *mols;
};
struct time_struct{
  int *id;
  double *time;
};


/* Parameters used to determine transition rates */
struct transitionParams{
  int array_size; //Equal to the number of levels (i.e NEQ)
  double *A_array; //Holds Einstein's As
  double *transition_rates; //Holds transition rates
  molData *md;
  int ispec;
  struct grid *gp;
  int *id;
  configInfo *par;
  struct time_struct time_struct;
  double *jbar_grid;
  int *nMaserWarnings;
};

/* Checks for errors when calling any CVode functions */
int check_retval(void *returnvalue, const char *funcname, int opt)
{
  int *retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
      funcname);
    return(1); }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int *) returnvalue;
    if (*retval < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with retval = %d\n\n",
        funcname, *retval);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
      funcname);
    return(1); }

  return(0);
}


/*ELECTRON TEMPERATURE FUNCTION*/
double Telec(double r, double Q, double Tkin){
  double Te,rcs;
  double Tmax = 1e4;
  rcs = 1.125e6 * pow(Q/1e29,0.75);
  if (r < rcs) {
    Te = Tkin;
  }
  else if (r > 2.*rcs) {
    Te = Tmax;
  }
  else {
    Te = Tkin + (Tmax - Tkin)*((r/rcs)-1.);
  }
  return Te;
}

/*ELECTRON DENSITY FUNCTION*/
double nelec(double r, double Q, double vexp, double Te, double rH, double xne){
  double Rrec, ne, krec, kion;
  kion = 4.1e-7;
  krec = 3e-13 * sqrt(300./Te); /* Recombination rate From RATE12, accountiong for cm3 to m3 conversion */
  Rrec = 3.2e6 * sqrt(Q/1e29);
  /* Equation 5 of Zakharov 2007 */
  ne = xne * sqrt(Q*kion/vexp/krec/(rH*rH)) * pow((Te/300.),0.15) * (Rrec/(r*r)) * (1-exp(-r/Rrec)) + (5e6/(rH*rH));
  
  return ne;
}



/*....................................................................*/
void
freeMolsWithBlends(struct molWithBlends *mols, const int numMolsWithBlends){
  int mi, li;

  if(mols != NULL){
    for(mi=0;mi<numMolsWithBlends;mi++){
      if(mols[mi].lines != NULL){
        for(li=0;li<mols[mi].numLinesWithBlends;li++)
          free(mols[mi].lines[li].blends);
        free(mols[mi].lines);
      }
    }
    free(mols);
  }
}

/*....................................................................*/
void
freeGridPointData(const int nSpecies, gridPointData *mol){
  /*
Note that this is called from within the multi-threaded block.
  */
  int i;
  if(mol!= NULL){
    for(i=0;i<nSpecies;i++){
      free(mol[i].jbar);
      free(mol[i].phot);
      free(mol[i].vfac);
      free(mol[i].vfac_loc);
    }
  }
}

/*....................................................................*/
void lineBlend(molData *m, configInfo *par, struct blendInfo *blends){
  /*
This obtains information on all the lines of all the radiating species which have other lines within some cutoff velocity separation.

A variable of type 'struct blendInfo' has a nested structure which can be illustrated diagrammaticaly as follows.

  Structs:  blendInfo   molWithBlends   lineWithBlends    blend

  Variables:  blends
      .numMolsWithBlends     ____________________
      .*mols--------------->|.molI               |
                            |.numLinesWithBlends |   ___________
                            |.*lines--------------->|.lineI     |
                            |____________________|  |.numBlends |           ________
                            |        etc         |  |.*blends------------->|.molJ   |
                                                    |___________|          |.lineJ  |
                                                    |    etc    |          |.deltaV |
                                                                           |________|
                                                                           |   etc  |

Pointers are indicated by a * before the attribute name and an arrow to the memory location pointed to.
  */
  int molI, lineI, molJ, lineJ;
  int nmwb, nlwb, numBlendsFound, li, bi;
  double deltaV;
  struct blend *tempBlends=NULL;
  struct lineWithBlends *tempLines=NULL;

  /* Dimension blends.mols first to the total number of species, then realloc later if need be.
  */
  (*blends).mols = malloc(sizeof(struct molWithBlends)*par->nSpecies);
  (*blends).numMolsWithBlends = 0;

  nmwb = 0;
  for(molI=0;molI<par->nSpecies;molI++){
    tempBlends = malloc(sizeof(struct blend)*m[molI].nline);
    tempLines  = malloc(sizeof(struct lineWithBlends)*m[molI].nline);

    nlwb = 0;
    for(lineI=0;lineI<m[molI].nline;lineI++){
      numBlendsFound = 0;
      for(molJ=0;molJ<par->nSpecies;molJ++){
        for(lineJ=0;lineJ<m[molJ].nline;lineJ++){
          if(!(molI==molJ && lineI==lineJ)){
            deltaV = (m[molJ].freq[lineJ] - m[molI].freq[lineI])*CLIGHT/m[molI].freq[lineI];
            if(fabs(deltaV)<maxBlendDeltaV){
              tempBlends[numBlendsFound].molJ   = molJ;
              tempBlends[numBlendsFound].lineJ  = lineJ;
              tempBlends[numBlendsFound].deltaV = deltaV;
              numBlendsFound++;
            }
          }
        }
      }

      if(numBlendsFound>0){
        tempLines[nlwb].lineI = lineI;
        tempLines[nlwb].numBlends = numBlendsFound;
        tempLines[nlwb].blends = malloc(sizeof(struct blend)*numBlendsFound);
        for(bi=0;bi<numBlendsFound;bi++)
          tempLines[nlwb].blends[bi] = tempBlends[bi];

        nlwb++;
      }
    }

    if(nlwb>0){
      (*blends).mols[nmwb].molI = molI;
      (*blends).mols[nmwb].numLinesWithBlends = nlwb;
      (*blends).mols[nmwb].lines = malloc(sizeof(struct lineWithBlends)*nlwb);
      for(li=0;li<nlwb;li++){
        (*blends).mols[nmwb].lines[li].lineI     = tempLines[li].lineI;
        (*blends).mols[nmwb].lines[li].numBlends = tempLines[li].numBlends;
        (*blends).mols[nmwb].lines[li].blends = malloc(sizeof(struct blend)*tempLines[li].numBlends);
        for(bi=0;bi<tempLines[li].numBlends;bi++)
          (*blends).mols[nmwb].lines[li].blends[bi] = tempLines[li].blends[bi];
      }

      nmwb++;
    }

    free(tempLines);
    free(tempBlends);
  }

  (*blends).numMolsWithBlends = nmwb;
  if(nmwb>0){
    if(!par->blend)
      if(!silent) warning("There are blended lines, but line blending is switched off.");

    (*blends).mols = realloc((*blends).mols, sizeof(struct molWithBlends)*nmwb);
  }else{
    if(par->blend)
      if(!silent) warning("Line blending is switched on, but no blended lines were found.");

    free((*blends).mols);
    (*blends).mols = NULL;
  }
}

/*....................................................................*/
void calcGridCollRates(configInfo *par, molData *md, struct grid *gp){
  int i,id,ipart,itrans,itemp,tnint=-1;
  struct cpData part;
  double fac;

  for(i=0;i<par->nSpecies;i++){
    for(id=0;id<par->ncell;id++){
      gp[id].mol[i].partner = malloc(sizeof(struct rates)*md[i].npart);
    }

    for(ipart=0;ipart<md[i].npart;ipart++){
      part = md[i].part[ipart];
      for(id=0;id<par->ncell;id++){
        for(itrans=0;itrans<part.ntrans;itrans++){
          if((gp[id].t[0]>part.temp[0])&&(gp[id].t[0]<part.temp[part.ntemp-1])){
            for(itemp=0;itemp<part.ntemp-1;itemp++){
              if((gp[id].t[0]>part.temp[itemp])&&(gp[id].t[0]<=part.temp[itemp+1])){
                tnint=itemp;
              }
            }
            fac=(gp[id].t[0]-part.temp[tnint])/(part.temp[tnint+1]-part.temp[tnint]);
            gp[id].mol[i].partner[ipart].t_binlow = tnint;
            gp[id].mol[i].partner[ipart].interp_coeff = fac;

    } else if(gp[id].t[0]<=part.temp[0]) {
      gp[id].mol[i].partner[ipart].t_binlow = 0;
      gp[id].mol[i].partner[ipart].interp_coeff = 0.0;
    } else {
      gp[id].mol[i].partner[ipart].t_binlow = part.ntemp-2;
      gp[id].mol[i].partner[ipart].interp_coeff = 1.0;
    }
        } /* End loop over transitions. */
      } /* End loop over grid points. */
    } /* End loop over collision partners. */
  } /* End loop over radiating molecules. */
}

/*....................................................................*/
void mallocGridCont(configInfo *par, molData *md, struct grid *gp){
  int id,si,li;

  for(id=0;id<par->ncell;id++){
    for(si=0;si<par->nSpecies;si++){
      gp[id].mol[si].cont = malloc(sizeof(*(gp[id].mol[si].cont))*md[si].nline);
      for(li=0;li<md[si].nline;li++){
        gp[id].mol[si].cont[li].dust = 0.0;
        gp[id].mol[si].cont[li].knu  = 0.0;
      }
    }
  }
}

/*....................................................................*/
void freeGridCont(configInfo *par, struct grid *gp){
  int id,si;

  for(id=0;id<par->ncell;id++){
    if(gp[id].mol==NULL)
      continue;

    for(si=0;si<par->nSpecies;si++){
      free(gp[id].mol[si].cont);
      gp[id].mol[si].cont = NULL;
    }
  }
}

/*....................................................................*/
void calcGridLinesDustOpacity(configInfo *par, molData *md, double *lamtab\
  , double *kaptab, const int nEntries, struct grid *gp){

  int iline,id,si;
  double *kappatab,gtd;
  gsl_spline *spline = NULL;
  gsl_interp_accel *acc = NULL;
  double *knus=NULL, *dusts=NULL;

  if(par->dust != NULL){
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline,nEntries);
    gsl_spline_init(spline,lamtab,kaptab,nEntries);
  }

  for(si=0;si<par->nSpecies;si++){
    kappatab = malloc(sizeof(*kappatab)*md[si].nline);
    knus     = malloc(sizeof(*knus)    *md[si].nline);
    dusts    = malloc(sizeof(*dusts)   *md[si].nline);

    if(par->dust == NULL){
      for(iline=0;iline<md[si].nline;iline++)
        kappatab[iline] = 0.;
    }else{
      for(iline=0;iline<md[si].nline;iline++)
        kappatab[iline] = interpolateKappa(md[si].freq[iline]\
                        , lamtab, kaptab, nEntries, spline, acc);
    }

    for(id=0;id<par->ncell;id++){
      gasIIdust(gp[id].x[0],gp[id].x[1],gp[id].x[2],&gtd);
      calcDustData(par, gp[id].dens, md[si].freq, gtd, kappatab, md[si].nline, gp[id].t, knus, dusts);
      for(iline=0;iline<md[si].nline;iline++){
        gp[id].mol[si].cont[iline].knu  = knus[iline];
        gp[id].mol[si].cont[iline].dust = dusts[iline];
      }
    }

    free(kappatab);
    free(knus);
    free(dusts);
  }

  if(par->dust != NULL){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
}

/*....................................................................*/
int
getNextEdge(double *inidir, const int startGi, const int presentGi\
  , struct grid *gp, const gsl_rng *ran){
  /*
The idea here is to select for the next grid point, that one which lies closest (with a little randomizing jitter) to the photon track, while requiring the direction of the edge to be in the 'forward' hemisphere of the photon direction.

Note that this is called from within the multi-threaded block.
  */
  int i,ni,niOfSmallest=-1,niOfNextSmallest=-1;
  double dirCos,distAlongTrack,dirFromStart[3],coord,distToTrackSquared,smallest=0.0,nextSmallest=0.0;
  const static double scatterReduction = 0.4;
  /*
This affects the ratio of N_2/N_1, where N_2 is the number of times the edge giving the 2nd-smallest distance from the photon track is chosen and N_1 ditto the smallest. Some ratio values obtained from various values of scatterReduction:

  scatterReduction  <N_2/N_1>
    1.0     0.42
    0.5     0.75
    0.4     0.90
    0.2     1.52

Note that the equivalent ratio value produced by the 1.6 code was 0.91.
  */

  i = 0;
  for(ni=0;ni<gp[presentGi].numNeigh;ni++){
    dirCos = dotProduct3D(inidir, gp[presentGi].dir[ni].xn);

    if(dirCos<=0.0)
  continue; /* because the edge points in the backward direction. */

    dirFromStart[0] = gp[presentGi].neigh[ni]->x[0] - gp[startGi].x[0];
    dirFromStart[1] = gp[presentGi].neigh[ni]->x[1] - gp[startGi].x[1];
    dirFromStart[2] = gp[presentGi].neigh[ni]->x[2] - gp[startGi].x[2];
    distAlongTrack = dotProduct3D(inidir, dirFromStart);

    coord = dirFromStart[0] - distAlongTrack*inidir[0];
    distToTrackSquared  = coord*coord;
    coord = dirFromStart[1] - distAlongTrack*inidir[1];
    distToTrackSquared += coord*coord;
    coord = dirFromStart[2] - distAlongTrack*inidir[2];
    distToTrackSquared += coord*coord;

    if(i==0){
      smallest = distToTrackSquared;
      niOfSmallest = ni;
    }else{
      if(distToTrackSquared<smallest){
        nextSmallest = smallest;
        niOfNextSmallest = niOfSmallest;
        smallest = distToTrackSquared;
        niOfSmallest = ni;
      }else if(i==1 || distToTrackSquared<nextSmallest){
        nextSmallest = distToTrackSquared;
        niOfNextSmallest = ni;
      }
    }

    i++;
  }

  /* Choose the edge to follow.
  */
  if(i>1){ /* then nextSmallest, niOfNextSmallest should exist. */
    if((smallest + scatterReduction*nextSmallest)*gsl_rng_uniform(ran)<smallest){
      return niOfNextSmallest;
    }else{
      return niOfSmallest;
    }
  }else if(i>0){
    return niOfSmallest;
  }else{
    if(!silent)
      bail_out("Photon propagation error - no valid edges.");
    exit(1);
  }
}

/*....................................................................*/
void calcLineAmpPWLin(struct grid *g, const int id, const int k\
  , const int molI, const double deltav, double *inidir, double *vfac_in, double *vfac_out){
  /*
Note that this is called from within the multi-threaded block.
  */

  /* convolution of a Gaussian with a box */
  double binv_this, binv_next, v[5];

  binv_this=g[id].mol[molI].binv;
  binv_next=(g[id].neigh[k])->mol[molI].binv;
  v[0]=deltav-dotProduct3D(inidir,g[id].vel);
  v[1]=deltav-dotProduct3D(inidir,&(g[id].v1[3*k]));
  v[2]=deltav-dotProduct3D(inidir,&(g[id].v2[3*k]));
  v[3]=deltav-dotProduct3D(inidir,&(g[id].v3[3*k]));
  v[4]=deltav-dotProduct3D(inidir,g[id].neigh[k]->vel);

  /* multiplying by the appropriate binv changes from velocity to doppler widths(?) */
  /* if the values were be no more than 2 erf table bins apart, we just take a single Gaussian */

  /*
  vfac_out is the lineshape for the part of the edge in the current Voronoi cell,
  vfac_in is for the part in the next cell
  */

  if (fabs(v[1]-v[0])*binv_this>(2.0*BIN_WIDTH)) {
     *vfac_out=0.5*geterf(v[0]*binv_this,v[1]*binv_this);
  } else *vfac_out=0.5*gaussline(0.5*(v[0]+v[1]),binv_this);
  if (fabs(v[2]-v[1])*binv_this>(2.0*BIN_WIDTH)) {
    *vfac_out+=0.5*geterf(v[1]*binv_this,v[2]*binv_this);
  } else *vfac_out+=0.5*gaussline(0.5*(v[1]+v[2]),binv_this);

  if (fabs(v[3]-v[2])*binv_next>(2.0*BIN_WIDTH)) {
     *vfac_in=0.5*geterf(v[2]*binv_next,v[3]*binv_next);
  } else *vfac_in=0.5*gaussline(0.5*(v[2]+v[3]),binv_next);
  if (fabs(v[4]-v[3])*binv_next>(2.0*BIN_WIDTH)) {
    *vfac_in+=0.5*geterf(v[3]*binv_next,v[4]*binv_next);
  } else *vfac_in+=0.5*gaussline(0.5*(v[3]+v[4]),binv_next);
}

/*....................................................................*/
void calcLineAmpLin(struct grid *g, const int id, const int k\
  , const int molI, const double deltav, double *inidir, double *vfac_in, double *vfac_out){
  /*
Note that this is called from within the multi-threaded block.
  */

  /* convolution of a Gaussian with a box */
  double binv_this, binv_next, v[3];

  binv_this=g[id].mol[molI].binv;
  binv_next=(g[id].neigh[k])->mol[molI].binv;
  v[0]=deltav-dotProduct3D(inidir,g[id].vel);
  v[2]=deltav-dotProduct3D(inidir,g[id].neigh[k]->vel);
  v[1]=0.5*(v[0]+v[2]);

  if (fabs(v[1]-v[0])*binv_this>(2.0*BIN_WIDTH)) {
     *vfac_out=geterf(v[0]*binv_this,v[1]*binv_this);
  } else *vfac_out+=gaussline(0.5*(v[0]+v[1]),binv_this);

  if (fabs(v[2]-v[1])*binv_next>(2.0*BIN_WIDTH)) {
     *vfac_in=geterf(v[1]*binv_next,v[2]*binv_next);
  } else *vfac_in+=gaussline(0.5*(v[1]+v[2]),binv_next);
}

/*....................................................................*/
void
calculateJBar(int id, struct grid *gp, molData *md, const gsl_rng *ran\
  , configInfo *par, const int nlinetot, struct blendInfo blends\
  , gridPointData *mp, double *halfFirstDs, int *nMaserWarnings){
  /*
Note that this is called from within the multi-threaded block.
  */

  int iphot,iline,here,there,firststep,neighI,numLinks=0;
  int nextMolWithBlend, nextLineWithBlend, molI, lineI, molJ, lineJ, bi;
  double segment,vblend_in,vblend_out,dtau,expDTau,ds_in=0.0,ds_out=0.0,pt_theta,pt_z,semiradius;
  double deltav[par->nSpecies],vfac_in[par->nSpecies],vfac_out[par->nSpecies],vfac_inprev[par->nSpecies];
  double expTau[nlinetot],inidir[3];
  double remnantSnu,velProj;
  char message[STR_LEN_0];

  for(iphot=0;iphot<gp[id].nphot;iphot++){
    firststep=1;
    iline = 0;
    for(molI=0;molI<par->nSpecies;molI++){
      for(lineI=0;lineI<md[molI].nline;lineI++){
        mp[molI].phot[lineI+iphot*md[molI].nline]=0.;
        expTau[iline]=1.;
        iline++;
      }
    }

    /* Choose random initial photon direction (the distribution used here is even over the surface of a sphere of radius 1).
    */
    pt_theta=gsl_rng_uniform(ran)*2*M_PI;
    pt_z=2*gsl_rng_uniform(ran)-1;
    semiradius = sqrt(1.-pt_z*pt_z);
    inidir[0]=semiradius*cos(pt_theta);
    inidir[1]=semiradius*sin(pt_theta);
    inidir[2]=pt_z;

    /* Choose the photon frequency/velocity offset.
    */
    segment=gsl_rng_uniform(ran)-0.5;
    /*
    Values of segment should be evenly distributed (considering the
    entire ensemble of photons) between -0.5 and +0.5.
    */

    for (molI=0;molI<par->nSpecies;molI++){
      /* Is factor 4.3=[-2.15,2.15] enough?? */
      deltav[molI]=4.3*segment*gp[id].mol[molI].dopb+dotProduct3D(inidir,gp[id].vel);
      /*
      This is the local (=evaluated at a grid point, not averaged over the local cell) lineshape.
      We store this for later use in ALI loops.
      */
      mp[molI].vfac_loc[iphot]=gaussline(deltav[molI]-dotProduct3D(inidir,gp[id].vel),gp[id].mol[molI].binv);
    }

    here = gp[id].id;

    /* Photon propagation loop */
    numLinks=0;
    while(!gp[here].sink){ /* Testing for sink at loop start is redundant for the first step, since we only start photons from non-sink points, but it makes for simpler code. */
      numLinks++;
      if(numLinks>par->ncell){
        if(!silent){
          snprintf(message, STR_LEN_0, "Bad grid? Too many links in photon path, point %d photon %d", id, iphot);
          bail_out(message);
        }
exit(1);
      }

      neighI = getNextEdge(inidir,id,here,gp,ran);

      there=gp[here].neigh[neighI]->id;

      if(firststep){
        firststep=0;
        ds_out=0.5*gp[here].ds[neighI]*dotProduct3D(inidir,gp[here].dir[neighI].xn);
        halfFirstDs[iphot]=ds_out;

        for(molI=0;molI<par->nSpecies;molI++){
          if(par->edgeVelsAvailable) {
            calcLineAmpPWLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
         } else
            calcLineAmpLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);

          mp[molI].vfac[iphot]=vfac_out[molI];
        }
        /*
        Contribution of the local cell to emission and absorption is done in updateJBar.
        We only store the vfac for the local cell for use in ALI loops.
        */
        here=there;
    continue;
      }

      /* If we've got to here, we have progressed beyond the first edge. Length of the new "in" edge is the length of the previous "out".
      */
      ds_in=ds_out;
      ds_out=0.5*gp[here].ds[neighI]*dotProduct3D(inidir,gp[here].dir[neighI].xn);

      for(molI=0;molI<par->nSpecies;molI++){
        vfac_inprev[molI]=vfac_in[molI];
        if(par->edgeVelsAvailable)
          calcLineAmpPWLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
        else
          calcLineAmpLin(gp,here,neighI,molI,deltav[molI],inidir,&vfac_in[molI],&vfac_out[molI]);
      }

      nextMolWithBlend = 0;
      iline = 0;
      for(molI=0;molI<par->nSpecies;molI++){
        nextLineWithBlend = 0;
        for(lineI=0;lineI<md[molI].nline;lineI++){
          double jnu_line_in=0., jnu_line_out=0., jnu_cont=0., jnu_blend=0.;
          double alpha_line_in=0., alpha_line_out=0., alpha_cont=0., alpha_blend=0.;

          sourceFunc_line(&md[molI],vfac_inprev[molI],&(gp[here].mol[molI]),lineI,&jnu_line_in,&alpha_line_in);
          sourceFunc_line(&md[molI],vfac_out[molI],&(gp[here].mol[molI]),lineI,&jnu_line_out,&alpha_line_out);
          sourceFunc_cont(gp[here].mol[molI].cont[lineI],&jnu_cont,&alpha_cont);

          /* cont and blend could use the same alpha and jnu counter, but maybe it's clearer this way */

          /* Line blending part.
          */
          if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI\
          && lineI==blends.mols[nextMolWithBlend].lines[nextLineWithBlend].lineI){

            for(bi=0;bi<blends.mols[nextMolWithBlend].lines[nextLineWithBlend].numBlends;bi++){
              molJ  = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].molJ;
              lineJ = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].lineJ;
              velProj = deltav[molI] - blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].deltaV;
        /*  */
              if(par->edgeVelsAvailable)
                calcLineAmpPWLin(gp,here,neighI,molJ,velProj,inidir,&vblend_in,&vblend_out);
              else
                calcLineAmpLin(gp,here,neighI,molJ,velProj,inidir,&vblend_in,&vblend_out);

        /* we should use also the previous vblend_in, but I don't feel like writing the necessary code now */
              sourceFunc_line(&md[molJ],vblend_out,&(gp[here].mol[molJ]),lineJ,&jnu_blend,&alpha_blend);
              /* note that sourceFunc* increment jnu and alpha, they don't overwrite it  */
            }

            nextLineWithBlend++;
            if(nextLineWithBlend>=blends.mols[nextMolWithBlend].numLinesWithBlends){
              nextLineWithBlend = 0;
              /* The reason for doing this is as follows. Firstly, we only enter the present IF block if molI has at least 1 line which is blended with others; and further, if we have now processed all blended lines for that molecule. Thus no matter what value lineI takes for the present molecule, it won't appear as blends.mols[nextMolWithBlend].lines[i].lineI for any i. Yet we will still test blends.mols[nextMolWithBlend].lines[nextLineWithBlend], thus we want nextLineWithBlend to at least have a sensible value between 0 and blends.mols[nextMolWithBlend].numLinesWithBlends-1. We could set nextLineWithBlend to any number in this range in safety, but zero is simplest. */
            }
          }
          /* End of line blending part */

    /* as said above, out-in split should be done also for blended lines... */

    dtau=(alpha_line_out+alpha_cont+alpha_blend)*ds_out;
          if(dtau < -MAX_NEG_OPT_DEPTH) dtau = -MAX_NEG_OPT_DEPTH;
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= (jnu_line_out+jnu_cont+jnu_blend)*ds_out;
          mp[molI].phot[lineI+iphot*md[molI].nline]+=expTau[iline]*remnantSnu;
    expTau[iline]*=expDTau;

    dtau=(alpha_line_in+alpha_cont+alpha_blend)*ds_in;
          if(dtau < -MAX_NEG_OPT_DEPTH) dtau = -MAX_NEG_OPT_DEPTH;
          calcSourceFn(dtau, par, &remnantSnu, &expDTau);
          remnantSnu *= (jnu_line_in+jnu_cont+jnu_blend)*ds_in;
          mp[molI].phot[lineI+iphot*md[molI].nline]+=expTau[iline]*remnantSnu;
    expTau[iline]*=expDTau;

          if(expTau[iline] > exp(MAX_NEG_OPT_DEPTH)){
            (*nMaserWarnings)++;
            expTau[iline]=exp(MAX_NEG_OPT_DEPTH);
          }

          iline++;
        } /* Next line this molecule. */

        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI)
          nextMolWithBlend++;
      }

      here=there;
    };

    /* Add cmb contribution.
    */
    iline = 0;
    for(molI=0;molI<par->nSpecies;molI++){
      for(lineI=0;lineI<md[molI].nline;lineI++){
        mp[molI].phot[lineI+iphot*md[molI].nline]+=expTau[iline]*md[molI].cmb[lineI];
        iline++;
      }
    }
  }
}

/*....................................................................*/
void
updateJBar(int posn, molData *md, struct grid *gp, const int molI\
  , configInfo *par, struct blendInfo blends, int nextMolWithBlend\
  , gridPointData *mp, double *halfFirstDs){
  /*
Note that this is called from within the multi-threaded block.
  */
  int lineI,iphot,bi,molJ,lineJ,nextLineWithBlend;
  double dtau,expDTau,remnantSnu,vsum=0.;
  
  for(lineI=0;lineI<md[molI].nline;lineI++) mp[molI].jbar[lineI]=0.;

  for(iphot=0;iphot<gp[posn].nphot;iphot++){
    if(mp[molI].vfac_loc[iphot]>0){
      nextLineWithBlend = 0;
      for(lineI=0;lineI<md[molI].nline;lineI++){
        double jnu=0.0;
        double alpha=0;

        sourceFunc_line(&md[molI],mp[molI].vfac[iphot],&(gp[posn].mol[molI]),lineI,&jnu,&alpha);
        sourceFunc_cont(gp[posn].mol[molI].cont[lineI],&jnu,&alpha);

        /* Line blending part.
        */
        if(par->blend && blends.mols!=NULL && molI==blends.mols[nextMolWithBlend].molI\
        && lineI==blends.mols[nextMolWithBlend].lines[nextLineWithBlend].lineI){
          for(bi=0;bi<blends.mols[nextMolWithBlend].lines[nextLineWithBlend].numBlends;bi++){
            molJ  = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].molJ;
            lineJ = blends.mols[nextMolWithBlend].lines[nextLineWithBlend].blends[bi].lineJ;
            /*
            The next line is not quite correct, because vfac may be different for other molecules due to different values of binv. Unfortunately we don't necessarily have vfac for molJ available yet.
            */
            sourceFunc_line(&md[molJ],mp[molI].vfac[iphot],&(gp[posn].mol[molJ]),lineJ,&jnu,&alpha);
      /* note that sourceFunc* increment jnu and alpha, they don't overwrite it  */
          }

          nextLineWithBlend++;
          if(nextLineWithBlend>=blends.mols[nextMolWithBlend].numLinesWithBlends){
            nextLineWithBlend = 0;
            /* The reason for doing this is as follows. Firstly, we only enter the present IF block if molI has at least 1 line which is blended with others; and further, if we have now processed all blended lines for that molecule. Thus no matter what value lineI takes for the present molecule, it won't appear as blends.mols[nextMolWithBlend].lines[i].lineI for any i. Yet we will still test blends.mols[nextMolWithBlend].lines[nextLineWithBlend], thus we want nextLineWithBlend to at least have a sensible value between 0 and blends.mols[nextMolWithBlend].numLinesWithBlends-1. We could set nextLineWithBlend to any number in this range in safety, but zero is simplest. */
          }
        }
        /* End of line blending part */

        dtau=alpha*halfFirstDs[iphot];
        calcSourceFn(dtau, par, &remnantSnu, &expDTau);
        remnantSnu *= jnu*halfFirstDs[iphot];
        mp[molI].jbar[lineI]+=mp[molI].vfac_loc[iphot]*(expDTau*mp[molI].phot[lineI+iphot*md[molI].nline]+remnantSnu);

      }
      vsum+=mp[molI].vfac_loc[iphot];
    }
  }
  for(lineI=0;lineI<md[molI].nline;lineI++) mp[molI].jbar[lineI] /= vsum;
}

/*....................................................................*/

void
getTransitionRates(molData *md, int ispec, struct grid *gp, int id, configInfo *par, int NEQ, double A[NEQ-1], double *p, double time, struct time_struct time_struct, double *jbar_grid, double *Pops, int *nMaserWarnings){
  int ipart,iline,k,l,ti, li, upper, lower, j;
  double radius, rnuc, Te, vexp, ne, aij, sigmaij, ve, bessel, ceij, gij, ceji, dens[md[ispec].npart];
  double jbar[par->pIntensity],molDens[par->nSpecies], tau, beta;

  rnuc = par->minScale;
  vexp = sqrt(gp[id].vel[0]* gp[id].vel[0]+gp[id].vel[1]*gp[id].vel[1]+gp[id].vel[2]*gp[id].vel[2]);
  radius = vexp*time + rnuc;
  density(radius,0.0,0.0,dens);

  if(time<time_struct.time[0])
    time = time_struct.time[0];
  else if(time>time_struct.time[par->pIntensity-1])
    time = time_struct.time[par->pIntensity-1];

  /* Initialize matrix with zeros */
  if(md[ispec].nlev<=0){
    if(!silent) bail_out("Matrix initialization error in solveStatEq");
    exit(1);
  }

  for(k=0;k<md[ispec].nlev;k++)
    for(l=0;l<md[ispec].nlev;l++)
      p[k * NEQ + l] = 0.0;

  /* Populate matrix with collisional transitions */
  for(ipart=0;ipart<md[ispec].npart;ipart++){
    struct cpData part = md[ispec].part[ipart];
    double *downrates = part.down;
    int di = md[ispec].part[ipart].densityIndex;
    if (di<0) continue;

    for(ti=0;ti<part.ntrans;ti++){
      int coeff_index = ti*part.ntemp + gp[id].mol[ispec].partner[ipart].t_binlow;
      double down = downrates[coeff_index]\
                  + gp[id].mol[ispec].partner[ipart].interp_coeff*(downrates[coeff_index+1]\
                  - downrates[coeff_index]);
      double up = down*md[ispec].gstat[part.lcu[ti]]/md[ispec].gstat[part.lcl[ti]]\
                *exp(-HCKB*(md[ispec].eterm[part.lcu[ti]]-md[ispec].eterm[part.lcl[ti]])/gp[id].t[0]);

      p[part.lcu[ti] * NEQ + part.lcl[ti]] = p[part.lcu[ti] * NEQ + part.lcl[ti]] + down*dens[ipart];
      p[part.lcl[ti] * NEQ + part.lcu[ti]] = p[part.lcl[ti] * NEQ + part.lcu[ti]] + up*dens[ipart];
    }

  }
  
  /*GENERATE ELECTRON COLLISIONAL RATES (only for gas 0) AND ADD TO MATRIX*/
  /*Formalism of Zakharov et al. (2007)*/
   Te = Telec(radius,par->Qwater,gp[id].t[0]);
   ne = nelec(radius,par->Qwater,vexp,Te,par->rHelio,par->xne);
   
   for(iline=0;iline<md[ispec].nline;iline++){
     aij = HPLANCK*md[ispec].freq[iline]/2./KBOLTZ/Te;
     sigmaij = ELEC_MASS*pow(ELEC_CHARGE,2)*pow(CLIGHT,3)*md[ispec].aeinst[iline]/16./pow(PI,2)/EPS_0/pow(HPLANCK,2)/pow(md[ispec].freq[iline],4);
     ve = sqrt(8.*KBOLTZ*Te/PI/ELEC_MASS);
     bessel = gsl_sf_bessel_K0(aij);
     ceij = ne*ve*sigmaij*2.*aij*exp(aij)*bessel;
     gij = md[ispec].gstat[md[ispec].lau[iline]]/md[ispec].gstat[md[ispec].lal[iline]];
     ceji = ne*ve*gij*sigmaij*2.*aij*exp(-aij)*bessel;

     p[md[ispec].lau[iline] * NEQ + md[ispec].lal[iline]] = p[md[ispec].lau[iline] * NEQ + md[ispec].lal[iline]] + ceij;
     p[md[ispec].lal[iline] * NEQ + md[ispec].lau[iline]] = p[md[ispec].lal[iline] * NEQ + md[ispec].lau[iline]] + ceji;

   }

   /* Add the pumping rates to the matrix */
   if(par->girdatfile!=NULL){
      for(k=0;k<md[ispec].nlev;k++){
        for(l=0;l<md[ispec].nlev;l++){
          if(k!=l){
            p[l * NEQ + k] = p[l * NEQ + k] + md[ispec].gir[l*md[ispec].nlev+k];
          }
        }
      }  
   }

 
   //Radiation trapping using the Escape Probaility method
   if(par->useEP){ 

   molNumDensity(radius,0.0,0.0, molDens); 
    for(li=0;li<md[ispec].nline;li++){
      upper=md[ispec].lau[li];
      lower=md[ispec].lal[li];

      //Calculating the optical depth
      tau = ((A[li]*pow(CLIGHT,3))/(8*PI*pow(md[ispec].freq[li],3))) * ((md[ispec].gstat[upper]/md[ispec].gstat[lower])*Pops[lower] - Pops[upper]) * ((molDens[ispec]* radius)/vexp);

      //If the optical depth is small, ignore it
      if (tau > -1.0e-6 && tau < 1.0e-6){
        beta = 1.0;
      }else if (tau>0.0) {
        beta = (2/(3*tau)) - exp(-tau/2)*(tau*(gsl_sf_bessel_Kn(2,tau/2)-gsl_sf_bessel_K1(tau/2))/3 -  gsl_sf_bessel_K1(tau/2)); 
      }else if (tau<0.0){ 
        beta = (1 - exp(-tau)) / tau;
      }    
        
      p[upper * NEQ + lower] = p[upper * NEQ + lower] + A[li]*beta;

    }//end for
  }else{
     
     for(li=0;li<md[ispec].nline;li++){
        upper=md[ispec].lau[li];
        lower=md[ispec].lal[li];
        p[upper * NEQ + lower] = p[upper * NEQ + lower] + A[li];
     }
  }
}

/*....................................................................*/
void lteOnePoint(molData *md, const int ispec, const double temp, double *pops){
  int ilev;
  double sum;

  sum = 0.0;
  for(ilev=0;ilev<md[ispec].nlev;ilev++){
    pops[ilev] = md[ispec].gstat[ilev]*exp(-HCKB*md[ispec].eterm[ilev]/temp);
    sum += pops[ilev];
  }
  for(ilev=0;ilev<md[ispec].nlev;ilev++)
    pops[ilev] /= sum;
}

/*....................................................................*/
void
LTE(configInfo *par, struct grid *gp, molData *md){
  int id,ispec;

  for(id=0;id<par->pIntensity;id++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      lteOnePoint(md, ispec, gp[id].t[0], gp[id].mol[ispec].pops);
    }
  }
  if(par->outputfile) popsout(par,gp,md);

}


/*....................................................................*/
/* Sets the Differential Equation (Pdot) to be solved by CVode */
int f(realtype t, N_Vector P, N_Vector Pdot, void *data){

  int i, j, NEQ, id;
  struct transitionParams *user_data = data;
  NEQ = user_data -> array_size;
  double *p = user_data -> transition_rates;
  id = *(user_data -> id);
  double Pops_array[NEQ];

  //Pasing P values to Pops_array for readability
  for(i=0; i < NEQ; ++i)
    Pops_array[i] = Ith(P,i);

  getTransitionRates(user_data->md,user_data->ispec,user_data->gp,id,user_data->par, NEQ, user_data -> A_array, p, t, user_data->time_struct, user_data->jbar_grid, Pops_array, user_data ->nMaserWarnings);

  //Initializing Pdot (otherwise it stores previous values between calls to CVode)
  // MAC: This seems unnecessary?
  for(i=0; i < NEQ; ++i)
    Ith(Pdot,i) = 0.0;

  //Setting Pdot
  for(i=0; i < NEQ; ++i)
    for(j=0; j < NEQ; ++j){
      if(i != j) 
        Ith(Pdot,i) = Ith(Pdot,i) + (Pops_array[j] * (p[j* NEQ +i])) - (Pops_array[i] * p[i * NEQ +j]);
    }
  return(0);
}

/*....................................................................*/
double getTime(struct grid *gp, int id, configInfo *par){
  double radius, time, rnuc, vexp;

  radius = sqrt(gp[id].x[0]* gp[id].x[0]+gp[id].x[1]*gp[id].x[1]+gp[id].x[2]*gp[id].x[2]);
  rnuc = par->minScale;
  vexp = sqrt(gp[id].vel[0]* gp[id].vel[0]+gp[id].vel[1]*gp[id].vel[1]+gp[id].vel[2]*gp[id].vel[2]);
  time = (radius - rnuc)/vexp;

  return(time);
}

/*....................................................................*/
void
solveStatEq(struct grid *gp, molData *md, const int ispec, configInfo *par\
  , struct blendInfo blends, int *nextMolWithBlend, gridPointData **mp\
  , double **halfFirstDs, int *nMaserWarnings){

  int id;
  realtype reltol, t;
  N_Vector P, abstol;
  SUNMatrix sunMatrix;
  SUNLinearSolver LS;
  void *cvode_mem;
  int i,j;
  int retval;

  gsl_vector *newpop = gsl_vector_alloc(md[ispec].nlev);

  /* Initializing parameters to be used by CVode */
  int NEQ = md[ispec].nlev; // Number of Equations .
  double p[NEQ][NEQ]; //Collisional rates
  double A[md[ispec].nline]; //Einstein As
  double Pops[NEQ]; //Level Populations
  double timearr[par->pIntensity], sorted_timearr[par->pIntensity]; //times through which the solver will iterate (each one corresponds to a specific gridpoint)

  for(id=0;id<md[ispec].nline;id++)
    A[id] = md[ispec].aeinst[id];
  
  for(id=0;id<par->pIntensity;id++){
    timearr[id] = getTime(gp,id,par);
    sorted_timearr[id] = timearr[id]; //Note: not sorted yet
  }
  //Sorting time array so that each grid point can be matched to its corresponding time
  qsort(sorted_timearr, par->pIntensity, sizeof(double), compare);

  struct time_struct time_struct;
  time_struct.id = malloc(sizeof(int)*par->pIntensity);
  time_struct.time = sorted_timearr;
  double current;

  for(i=0;i<par->pIntensity;i++){
    current = time_struct.time[i];
    for(j=0;j<par->pIntensity;j++) //TODO: More efficient algorithm than sequential search could be implemented
      if(current==timearr[j])
        time_struct.id[i] = j; //holds sorted ids according to the time(radius) of its corresponding gridpoint
  }
  double (*jbar_grid)[md[ispec].nline] = malloc(sizeof(double[par->pIntensity][md[ispec].nline]));

  /*Initializing Pops */
  for(i=0;i<md[ispec].nlev;i++)
    Pops[i] = gp[time_struct.id[0]].mol[ispec].pops[i]; //we use gp[time_struct.id[0]] since we only need to initialize the level populations for the initial time

  struct transitionParams user_data = {NEQ, A, *p, md,ispec,gp,time_struct.id,par,time_struct, *jbar_grid, nMaserWarnings}; 

  P = abstol = NULL;
  sunMatrix = NULL;
  LS = NULL;
  cvode_mem = NULL;

  /* Create serial vector of length NEQ for I.C. and abstol */
  P = N_VNew_Serial(NEQ);
  if (check_retval((void *)P, "N_VNew_Serial", 0)) return;
  abstol = N_VNew_Serial(NEQ); 
  if (check_retval((void *)abstol, "N_VNew_Serial", 0)) return;

  /* Set the scalar relative tolerance */
  reltol = RTOL;

  /* Set the initial population values and the tolerances*/
  for(i=0; i < NEQ; ++i){
    Ith(P,i) = Pops[i];
    Ith(abstol,i) = ATOL;
  }

  /* Call CVodeCreate to create the solver memory and specify the 
  * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF);
  if (check_retval((void *)cvode_mem, "CVodeCreate", 0)) return;

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f,  time_struct.time[0], P);
  if (check_retval(&retval, "CVodeInit", 1)) return;

  retval = CVodeSetUserData(cvode_mem, &user_data);
  if (check_retval(&retval, "CVodeSetUserData", 0)) return;

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return;

  retval = CVodeSetMaxNumSteps(cvode_mem, 5000);
    if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) return;

  /* Create dense SUNMatrix for use in linear solves */
  sunMatrix = SUNDenseMatrix(NEQ, NEQ);
  if(check_retval((void *)sunMatrix, "SUNDenseMatrix", 0)) return;

  /* Create dense SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_Dense(P, sunMatrix);
  if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return;

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, sunMatrix);
  if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return;

  /* Set the user-supplied Jacobian routine Jac */
  //retval = CVodeSetJacFn(cvode_mem, Jac);
  //if(check_retval(&retval, "CVodeSetJacFn", 1)) return;

  /* Use a difference quotient Jacobian */
  retval = CVodeSetJacFn(cvode_mem, NULL);
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return;

  //We start at i=1 (instead of i=0) because time_struct.time[i] indicates the first output time, and we are starting the run from time = time_struct.time[0]
  for (i = 1; i < par->pIntensity; i++){
    retval = 1;
    fflush(stdout);
    while(retval!=0){
      retval = CVode(cvode_mem, time_struct.time[i], P, &t, CV_NORMAL);
      if(retval==-3) printf("Continuing anyway (check populations!)\n");
    }
    
    if(retval == CV_SUCCESS){
      for(j=0; j < NEQ; ++j) 
        Pops[j] = Ith(P,j);
        

      for(j=0;j<md[ispec].nlev;j++){ 
        gsl_vector_set(newpop,j,Pops[j]);
        gsl_vector_set(newpop,j,gsl_max(gsl_vector_get(newpop,j),EPS)); 
        gp[time_struct.id[i]].mol[ispec].pops[j]= gsl_vector_get(newpop,j); 
      }
      user_data.id++;
    } 
   
  }
  
  /* Free P and abstol vectors */
  N_VDestroy(P);
  N_VDestroy(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  /* Free the linear solver memory */
  SUNLinSolFree(LS);

  /* Free the matrix memory */
  SUNMatDestroy(sunMatrix);
  
  gsl_vector_free(newpop);
  free(jbar_grid);
}


/*....................................................................*/
int
levelPops(molData *md, configInfo *par, struct grid *gp, int *popsdone, double *lamtab, double *kaptab, const int nEntries){

  int id,ispec,i,nVerticesDone,nItersDone,nlinetot;
  int totalNMaserWarnings=0;
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;
  struct blendInfo blends;
  char message[STR_LEN_0];
  int RNG_seeds[par->nThreads];
  gsl_error_handler_t *defaultErrorHandler=NULL;
  int nextMolWithBlend[par->pIntensity],nMaserWarnings[par->pIntensity];
  for(id=0;id<par->pIntensity;id++){
    nMaserWarnings[id] = 0;
  }

  nlinetot = 0;
  for(ispec=0;ispec<par->nSpecies;ispec++)
    nlinetot += md[ispec].nline;

  if(par->lte_only){
    LTE(par,gp,md);
    if(par->outputfile) popsout(par,gp,md);
  }

   /* Non-LTE */
  else{

    /* Random number generator */
    gsl_rng *ran = gsl_rng_alloc(ranNumGenType);
    if(fixRandomSeeds)
      gsl_rng_set(ran, 1237106) ;
    else 
      gsl_rng_set(ran,time(0));

    gsl_rng **threadRans;
    threadRans = malloc(sizeof(gsl_rng *)*par->pIntensity);

    for (i=0;i<par->pIntensity;i++){
      threadRans[i] = gsl_rng_alloc(ranNumGenType);
      if (par->resetRNG==1) RNG_seeds[i] = (int)(gsl_rng_uniform(ran)*1e6);
      else gsl_rng_set(threadRans[i],(int)(gsl_rng_uniform(ran)*1e6));
    }

   calcGridCollRates(par,md,gp);
   freeGridCont(par, gp);
   mallocGridCont(par, md, gp);
   calcGridLinesDustOpacity(par, md, lamtab, kaptab, nEntries, gp);

   /* Check for blended lines */
   lineBlend(md, par, &blends);

   /* Initialize populations with Boltzmann distribution (assuming LTE) */
   LTE(par,gp,md);

   if(par->outputfile) popsout(par,gp,md);

   defaultErrorHandler = gsl_set_error_handler_off();


  gridPointData *mp[par->pIntensity];
  double *halfFirstDs[par->pIntensity];

  calcGridMolSpecNumDens(par,md,gp);
  totalNMaserWarnings = 0;
  nVerticesDone=0;

  //TODO: This for loop could be parallelized
  for(id=0;id<par->pIntensity;id++){
    ++nVerticesDone;
    nMaserWarnings[id]=0;
    nextMolWithBlend[id] = 0;
    mp[id]=malloc(sizeof(gridPointData)*par->nSpecies);
    halfFirstDs[id] = malloc(sizeof(*halfFirstDs)*gp[id].nphot);

    for (ispec=0;ispec<par->nSpecies;ispec++){
      mp[id][ispec].jbar = malloc(sizeof(double)*md[ispec].nline);
      mp[id][ispec].phot = malloc(sizeof(double)*md[ispec].nline*gp[id].nphot);
      mp[id][ispec].vfac = malloc(sizeof(double)*                gp[id].nphot);
      mp[id][ispec].vfac_loc = malloc(sizeof(double)*            gp[id].nphot);
    }
    if(gp[id].dens[0] < 0 && gp[id].t[0] < 0){
      printf("\nError on grid point = %d\n, aborting", id);
      exit(1);
    }

  }
  for(ispec=0;ispec<par->nSpecies;ispec++){
    solveStatEq(gp,md,ispec,par,blends,nextMolWithBlend,mp,halfFirstDs, nMaserWarnings); 
    for(i=0;i<par->pIntensity;i++)
      if(par->blend && blends.mols!=NULL && ispec==blends.mols[nextMolWithBlend[i]].molI)
        nextMolWithBlend[i] = nextMolWithBlend[i] + 1;
  }

  for(id=0;id<par->pIntensity;id++){
    totalNMaserWarnings = nMaserWarnings[id];
  }

  if(!silent && totalNMaserWarnings>0){
    snprintf(message, STR_LEN_0, "Maser warning: optical depth dropped below -%4.1f %d times this iteration.", MAX_NEG_OPT_DEPTH, totalNMaserWarnings);
    warning(message);
  }

  if(!silent) warning("");
  for (i=0;i<par->pIntensity;i++){
    freeGridPointData(par->nSpecies, mp[i]);
    free(halfFirstDs[i]);
  }
  if(par->outputfile != NULL) popsout(par,gp,md);

  freeMolsWithBlends(blends.mols, blends.numMolsWithBlends);
  freeGridCont(par, gp);
  for (i=0;i<par->pIntensity;i++)
    gsl_rng_free(threadRans[i]);
  free(threadRans);
  gsl_rng_free(ran);
 }//end else

  par->dataFlags |= (1 << DS_bit_populations);

  if(par->binoutputfile != NULL) binpopsout(par,gp,md);

  *popsdone=1;

  return (1);
}