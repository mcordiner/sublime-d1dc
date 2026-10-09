/*
 *  solver.c
 *  
 *  This file is part of SUBLIMED1D (SUBlimating gases in LIME, version D1D)
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *  Copyright (C) 2023 Martin Cordiner, Emmanuel Garcia-Berrios and Kristen Darnell (NASA GSFC)
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver      */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */

#include "sublime.h"
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
struct cvode_physdata{
  int array_size; //Equal to the number of levels (i.e NEQ)
  double *A_array; //Holds Einstein's As
  molData *md;
  int ispec;
  struct grid *gp;
  configInfo *par;
  double *jbar_grid;
  int *nMaserWarnings;
  struct CKCdata CKCdata;

  /* --- Caching & heap-allocated working memory (ported from 3D solver) ---
     p_rates/coll_rates were previously a stack-based "double p[NEQ][NEQ]"
     VLA. Since NEQ scales with the number of energy levels, doubling the
     line list quadruples that array's footprint, which was blowing the
     stack and causing the segfault. Both are now malloc'd once per call to
     solveStatEq() and reused across every CVode call to f(). */
  double *p_rates;      // Working rate matrix, rebuilt (cheaply) every f() call
  double *coll_rates;   // Radius-dependent collisional/electron/pump rates, cached
  double last_radius;   // Radius at which coll_rates was last computed
  double vexp;          // Cached expansion velocity at last_radius
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

/* ELECTRON TEMPERATURE FUNCTION: ORIGINAL VERSION */
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

/* ELECTRON DENSITY FUNCTION: ORIGINAL VERSION */
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
Note that this is called from within the multi-threaded block.
  */
  int i,ni,niOfSmallest=-1,niOfNextSmallest=-1;
  double dirCos,distAlongTrack,dirFromStart[3],coord,distToTrackSquared,smallest=0.0,nextSmallest=0.0;
  const static double scatterReduction = 0.4;

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

  if(i>1){
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
  double binv_this, binv_next, v[5];

  binv_this=g[id].mol[molI].binv;
  binv_next=(g[id].neigh[k])->mol[molI].binv;
  v[0]=deltav-dotProduct3D(inidir,g[id].vel);
  v[1]=deltav-dotProduct3D(inidir,&(g[id].v1[3*k]));
  v[2]=deltav-dotProduct3D(inidir,&(g[id].v2[3*k]));
  v[3]=deltav-dotProduct3D(inidir,&(g[id].v3[3*k]));
  v[4]=deltav-dotProduct3D(inidir,g[id].neigh[k]->vel);

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
getTransitionRates(struct cvode_physdata *user_data, double radius, double *Pops){

  int itemp,ipart,t_binlow,iline,k,l,ti,li,upper,lower,tnint=-1;
  double rnuc, Te, ne, aij, sigmaij, ve, bessel, ceij, gij, ceji;
  double interp_coeff, vkin, collRate, tau, beta;

  molData *md      = user_data->md;
  int ispec        = user_data->ispec;
  configInfo *par  = user_data->par;
  struct grid *gp  = user_data->gp;
  int NEQ          = user_data->array_size;
  double *A        = user_data->A_array;
  double *p        = user_data->p_rates;
  double *coll_rates = user_data->coll_rates;

  /* 1. ONLY RECALCULATE COLLISIONAL/ELECTRON/PUMP RATES IF RADIUS HAS CHANGED
        (Optimization ported from the 3D solver). CVode frequently calls f()
        - and hence this routine - multiple times at the same radius (e.g.
        during Newton iterations / finite-difference Jacobian evaluation), so
        this cache avoids a lot of redundant O(NEQ^2) work. This becomes more
        valuable the larger NEQ (i.e. nlev) gets. */
  if (radius != user_data->last_radius) {

    double dens[md[ispec].npart], tkin[md[ispec].npart];
    double LTEpops[md[ispec].nlev];
    double vel[DIM];

    rnuc = par->minScale;
    density(par,radius,0.0,0.0,dens);
    temperature(par,radius,0.0,0.0,tkin);
    lteOnePoint(md, ispec, tkin[0], LTEpops);

    velocity(par,0.,0.,radius,vel);
    user_data->vexp = sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);

    /* Initialize cached matrix with zeros */
    if(md[ispec].nlev<=0){
      if(!silent) bail_out("Matrix initialization error in solveStatEq");
      exit(1);
    }
    for(k=0;k<NEQ*NEQ;k++) coll_rates[k] = 0.0;

    /* Populate matrix with collisional transitions */
    for(ipart=0;ipart<md[ispec].npart;ipart++){
      struct cpData part = md[ispec].part[ipart];
      double *downrates = part.down;
      int di = md[ispec].part[ipart].densityIndex;
      if (di<0) continue;

      /* Collision temperature interpolation coefficients */
      if((tkin[ipart]>part.temp[0])&&(tkin[ipart]<part.temp[part.ntemp-1])){
        for(itemp=0;itemp<part.ntemp-1;itemp++){
          if((tkin[ipart]>part.temp[itemp])&&(tkin[ipart]<=part.temp[itemp+1])){
            tnint=itemp;
          }
        }
        interp_coeff =(tkin[ipart]-part.temp[tnint])/(part.temp[tnint+1]-part.temp[tnint]);
        t_binlow = tnint;

      } else if(tkin[ipart]<=part.temp[0]) {
        t_binlow = 0;
        interp_coeff = 0.0;
      } else {
        t_binlow = part.ntemp-2;
        interp_coeff = 1.0;
      }

      if (par->xsec == 0){
        for(ti=0;ti<part.ntrans;ti++){
          int coeff_index = ti*part.ntemp + t_binlow;
          double down = downrates[coeff_index]
                      + interp_coeff*(downrates[coeff_index+1] - downrates[coeff_index]);
          double up = down*md[ispec].gstat[part.lcu[ti]]/md[ispec].gstat[part.lcl[ti]]
                    *exp(-HCKB*(md[ispec].eterm[part.lcu[ti]]-md[ispec].eterm[part.lcl[ti]])/tkin[ipart]);

          coll_rates[part.lcu[ti] * NEQ + part.lcl[ti]] += down*dens[ipart];
          coll_rates[part.lcl[ti] * NEQ + part.lcu[ti]] += up*dens[ipart];
        }
      }else{
        /* If xsec was provided, then use the Meudon approximation */
        vkin = sqrt(8.0*KBOLTZ*tkin[ipart]/PI * (1.0/md[ispec].amass + 1.0/MATM));
        collRate = dens[ipart]*vkin*par->xsec;
        /* CACHE OPTIMIZATION: loop order swapped (l outer, k inner) so that
           coll_rates[l*NEQ+k] is accessed contiguously along the inner loop,
           matching row-major storage - avoids cache misses, esp. important
           as NEQ (nlev) grows large. */
        for(l=0;l<md[ispec].nlev;l++){
          for(k=0;k<md[ispec].nlev;k++){
            coll_rates[l * NEQ + k] += collRate * LTEpops[k];
          }
        }
      }
    }

    /*GENERATE ELECTRON COLLISIONAL RATES AND ADD TO MATRIX*/
    //Presently, only electrons produced from collision parter 0 are considered.
    //It would be easy enough to add others, but the partner production rates would be needed as an input
    //parameter, like par->Qpartner, and their temperatures can be given as tkin[n] from temperature()

    if (par->useCKCdata == 1 || par->useCKCdata == 2) {
      /* Interpolating from the CKC code */
      Te = get_Telec (&user_data->CKCdata, par->Qwater, par->rHelio, radius);
      ne = get_nelec (&user_data->CKCdata, par->Qwater, par->rHelio, radius);
      //Below is for debugging only
      FILE *fPtr;
      char *filePath = "output/elec_data.txt";
      fPtr = fopen(filePath, "a");
      fprintf(fPtr, "%12.3e %12.3e %12.3e\n", radius, Te, ne);
      fclose(fPtr);
    }
    else {
      /*Formalism of Zakharov et al. (2007)*/
      Te = Telec(radius,par->Qwater,tkin[0]);
      ne = nelec(radius,par->Qwater,user_data->vexp,Te,par->rHelio,par->xne);
    }

    for(iline=0;iline<md[ispec].nline;iline++){
      aij = HPLANCK*md[ispec].freq[iline]/2./KBOLTZ/Te;
      sigmaij = ELEC_MASS*pow(ELEC_CHARGE,2)*pow(CLIGHT,3)*md[ispec].aeinst[iline]/16./pow(PI,2)/EPS_0/pow(HPLANCK,2)/pow(md[ispec].freq[iline],4);
      ve = sqrt(8.*KBOLTZ*Te/PI/ELEC_MASS);
      bessel = gsl_sf_bessel_K0(aij);
      ceij = ne*ve*sigmaij*2.*aij*exp(aij)*bessel;
      gij = md[ispec].gstat[md[ispec].lau[iline]]/md[ispec].gstat[md[ispec].lal[iline]];
      ceji = ne*ve*gij*sigmaij*2.*aij*exp(-aij)*bessel;

      coll_rates[md[ispec].lau[iline] * NEQ + md[ispec].lal[iline]] += ceij;
      coll_rates[md[ispec].lal[iline] * NEQ + md[ispec].lau[iline]] += ceji;
    }

    /* Add the pumping rates to the matrix */
    if(par->girdatfile!=NULL){
      /* CACHE OPTIMIZATION: l outer, k inner, matching row-major storage */
      for(l=0;l<md[ispec].nlev;l++){
        for(k=0;k<md[ispec].nlev;k++){
          if(k!=l){
            coll_rates[l * NEQ + k] += md[ispec].gir[l*md[ispec].nlev+k];
          }
        }
      }
    }

    user_data->last_radius = radius;
  } /* end of "if radius changed" cached block */

  /* 2. COPY CACHED RATES INTO THE WORKING MATRIX. This copy is still O(NEQ^2)
        every call, but it's much cheaper than recomputing collision rates,
        electron rates and interpolation coefficients from scratch, and it
        correctly allows per-call terms (escape probability / optical depth,
        which depend on the current Pops vector, not just radius) to be
        added fresh below without corrupting the cached base rates. */
  for(k=0;k<NEQ*NEQ;k++) p[k] = coll_rates[k];

  //Radiation trapping using the Escape Probability method
  if(par->useEP){

    double molDens[par->nSpecies];
    molNumDensity(par, radius,0.0,0.0, molDens);

    for(li=0;li<md[ispec].nline;li++){
      upper=md[ispec].lau[li];
      lower=md[ispec].lal[li];

      //Calculating the optical depth
      tau = ((A[li]*pow(CLIGHT,3))/(8*PI*pow(md[ispec].freq[li],3))) * ((md[ispec].gstat[upper]/md[ispec].gstat[lower])*Pops[lower] - Pops[upper]) * ((molDens[ispec]* radius)/user_data->vexp);

      //If the optical depth is small, ignore it
      if (tau > -1.0e-6 && tau < 1.0e-6){
        beta = 1.0;
      }else if (tau>0.0) {
        beta = (2/(3*tau)) - exp(-tau/2)*(tau*(gsl_sf_bessel_Kn(2,tau/2)-gsl_sf_bessel_K1(tau/2))/3 -  gsl_sf_bessel_K1(tau/2));
      }else if (tau<0.0){
        beta = (1 - exp(-tau)) / tau;
      }

      p[upper * NEQ + lower] += A[li]*beta;

    }//end for
  }else{

    for(li=0;li<md[ispec].nline;li++){
      upper=md[ispec].lau[li];
      lower=md[ispec].lal[li];
      p[upper * NEQ + lower] += A[li];
    }
  }
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
}


/*....................................................................*/
/* Sets the Differential Equation (Pdot) to be solved by CVode */
int f(realtype radius, N_Vector P, N_Vector Pdot, void *data){

  int i, j, NEQ;
  struct cvode_physdata *user_data = data;
  NEQ = user_data -> array_size;
  double *p = user_data -> p_rates;
  double Pops_array[NEQ];

  /* Pasing P values to Pops_array for readability */
  for(i=0; i < NEQ; ++i)
    Pops_array[i] = Ith(P,i);

  getTransitionRates(user_data, radius, Pops_array);

  /* Initializing Pdot (otherwise it stores previous values between calls to CVode) */
  for(i=0; i < NEQ; ++i)
    Ith(Pdot,i) = 0.0;

  /* Cache-friendly matrix-vector product, split into two passes so that
     both inner loops walk p[] contiguously in row-major order rather than
     striding across rows - this was previously a single combined loop
     indexing p[j*NEQ+i] and p[i*NEQ+j] together, which is cache-unfriendly
     and gets progressively worse as NEQ (nlev) grows. */

  /* 1. Add population flowing FROM j TO i */
  for(j=0; j < NEQ; ++j){
    double pop_j = Pops_array[j];
    for(i=0; i < NEQ; ++i){
      if(i != j){
        Ith(Pdot,i) += pop_j * p[j * NEQ + i];
      }
    }
  }

  /* 2. Subtract population flowing FROM i TO j */
  for(i=0; i < NEQ; ++i){
    double pop_i = Pops_array[i];
    double out_rate_sum = 0.0;
    for(j=0; j < NEQ; ++j){
      if(i != j){
        out_rate_sum += p[i * NEQ + j];
      }
    }
    Ith(Pdot,i) -= pop_i * out_rate_sum;
  }

  /* Multiply equations by 1/V to get dP/dr */
  double inv_vexp = 1.0 / user_data->vexp;
  for(i=0; i < NEQ; ++i){
    Ith(Pdot,i) *= inv_vexp;
  }

  return(0);
}

void reduceTol(N_Vector abstol, realtype *reltol, void *cvode_mem, double factor, int NEQ){
   printf("INFO: Reducing RTOL and ATOL by 0.1 and restarting CVODE\n");
   int i, retval;

   for(i=0; i < NEQ; ++i){
      Ith(abstol,i) = Ith(abstol,i) * factor;
   }

   *reltol = (*reltol) * factor;

   retval = CVodeSVtolerances(cvode_mem, *reltol, abstol);
   if (check_retval(&retval, "CVodeSVtolerances", 1)) return;
}

/*....................................................................*/
/* Exact Jacobian-times-vector product, for use with CVodeSetJacTimes.
   VALID ONLY WHEN par->useEP == 0. In that case the rate
   matrix p is independent of the population vector P (f(t,P) is exactly
   linear: f = (1/vexp)*p^T*P), so J = (1/vexp)*p^T
*/

int JacTimesVec(N_Vector v, N_Vector Jv, realtype t, N_Vector P, N_Vector fP,
                 void *user_data_ptr, N_Vector tmp){

  struct cvode_physdata *user_data = user_data_ptr;
  int NEQ   = user_data->array_size;
  double *p = user_data->p_rates;
  int i,j;

  for(i=0; i < NEQ; ++i)
    Ith(Jv,i) = 0.0;

  /* Contribution flowing FROM j TO i */
  for(j=0; j < NEQ; ++j){
    double v_j = Ith(v,j);
    for(i=0; i < NEQ; ++i){
      if(i != j)
        Ith(Jv,i) += v_j * p[j * NEQ + i];
    }
  }

  /* Contribution flowing FROM i TO j (subtracted) */
  for(i=0; i < NEQ; ++i){
    double out_rate_sum = 0.0;
    for(j=0; j < NEQ; ++j){
      if(i != j)
        out_rate_sum += p[i * NEQ + j];
    }
    Ith(Jv,i) -= Ith(v,i) * out_rate_sum;
  }

  double inv_vexp = 1.0 / user_data->vexp;
  for(i=0; i < NEQ; ++i)
    Ith(Jv,i) *= inv_vexp;

  return 0;
}

/*....................................................................*/
void
solveStatEq(struct grid *gp, molData *md, const int ispec, configInfo *par\
  , struct blendInfo blends, int *nextMolWithBlend, gridPointData **mp\
  , double **halfFirstDs, int *nMaserWarnings, double *radii, int *gp_sorter){

  int id;
  realtype reltol, t;
  N_Vector P, abstol;
  SUNMatrix sunMatrix = NULL;
  SUNLinearSolver LS = NULL;
  void *cvode_mem;
  int i,j, cvodeErrs = 0;
  int retval,cvstatus;
  int using_dense; /* Tier tracker for solver fallback */

  /* Initializing parameters to be used by CVode */
  int NEQ = md[ispec].nlev; // Number of Equations
  double A[md[ispec].nline]; //Einstein As
  double Pops[NEQ]; //Level Populations

  for(id=0;id<md[ispec].nline;id++)
    A[id] = md[ispec].aeinst[id];

  double (*jbar_grid)[md[ispec].nline] = malloc(sizeof(double[par->pIntensity][md[ispec].nline]));

  /*Initializing Pops */
  for(i=0;i<md[ispec].nlev;i++)
    Pops[i] = gp[gp_sorter[0]].mol[ispec].pops[i];

  struct CKCdata CKCdata;
  if (par->useCKCdata == 1) {
    readCKCdata(&CKCdata);
  }
  else if (par->useCKCdata == 2) {
    readCKCfile(&CKCdata, par->CKCTeFile, par->CKCneFile, par->Qwater, par->rHelio);
  }

  /* --- HEAP ALLOCATION FIX ---
     p_rates and coll_rates were previously "double p[NEQ][NEQ]", a stack
     VLA sized by nlev^2. With a doubled line list this quadruples in size
     and was overflowing the stack almost immediately on entry to this
     function (confirmed by the EXC_BAD_ACCESS/code=2 crash in solveStatEq,
     consistent with hitting the stack guard page). Both are now malloc'd
     once per call and reused across the entire CVODE integration, matching
     the approach already used successfully in the 3D solver. */
  double *p_rates    = malloc(sizeof(double) * NEQ * NEQ);
  double *coll_rates = malloc(sizeof(double) * NEQ * NEQ);
  if (check_retval((void *)p_rates, "malloc p_rates", 2)) return;
  if (check_retval((void *)coll_rates, "malloc coll_rates", 2)) return;

  struct cvode_physdata user_data = {
      NEQ, A, md, ispec, gp, par, *jbar_grid, nMaserWarnings, CKCdata,
      p_rates, coll_rates, -1.0, 1.0 /* last_radius, vexp initializers */
  };

  P = abstol = NULL;
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

  retval = CVodeInit(cvode_mem, f,  radii[0], P);
  if (check_retval(&retval, "CVodeInit", 1)) return;

  retval = CVodeSetUserData(cvode_mem, &user_data);
  if (check_retval(&retval, "CVodeSetUserData", 0)) return;

  retval = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_retval(&retval, "CVodeSVtolerances", 1)) return;

  retval = CVodeSetMaxNumSteps(cvode_mem, 5000);
  if (check_retval(&retval, "CVodeSetMaxNumSteps", 1)) return;

/* --- SOLVER SELECTION ---
     Force Dense whenever useEP>0: the escape-probability term makes the
     rate matrix p depend on Pops (nonlinear system), which invalidates
     JacTimesVec above, and we have not validated SPGMR's convergence
     behaviour for that nonlinear case. For useEP==0 (confirmed linear
     system), prefer SPGMR above DENSE_NEQ_CUTOFF, backed by the exact
     Jacobian-vector product and a widened Krylov subspace (20, rather
     than SUNDIALS' small default) to guard against the convergence
     failures that made SPGMR's earlier (FD-Jacobian) behaviour unreliable. */
  int forceDense = (par->useEP > 0);

  if (!forceDense && NEQ > DENSE_NEQ_CUTOFF) {
    using_dense = 0;
    LS = SUNLinSol_SPGMR(P, PREC_NONE, 20);
    if(check_retval((void *)LS, "SUNLinSol_SPGMR", 0)) return;

    retval = CVodeSetLinearSolver(cvode_mem, LS, NULL);
    if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return;

    retval = CVodeSetJacTimes(cvode_mem, NULL, JacTimesVec);
    if(check_retval(&retval, "CVodeSetJacTimes", 1)) return;
  } else {
    using_dense = 1;
    sunMatrix = SUNDenseMatrix(NEQ, NEQ);
    if(check_retval((void *)sunMatrix, "SUNDenseMatrix", 0)) return;

    LS = SUNLinSol_Dense(P, sunMatrix);
    if(check_retval((void *)LS, "SUNLinSol_Dense", 0)) return;

    retval = CVodeSetLinearSolver(cvode_mem, LS, sunMatrix);
    if(check_retval(&retval, "CVodeSetLinearSolver", 1)) return;
  }

  printf("Linear solver selected for species %d: %s (NEQ=%d, useEP=%d)\n",
         ispec, using_dense ? "Dense" : "Krylov", NEQ, par->useEP);
  fflush(stdout);

  retval = CVodeSetJacFn(cvode_mem, NULL);
  if(check_retval(&retval, "CVodeSetJacFn", 1)) return;

  printf("Starting CVODE time-dependent solver...\n");
  fflush(stdout);

  // Call CVODE for each radius
  // Do/while loop to catch CVODE error test failure status
  do{
     for(i=1; i<par->pIntensity; i++){
         cvstatus = CVode(cvode_mem, radii[i], P, &t, CV_NORMAL);

       if(cvstatus == CV_SUCCESS){
         for(j=0;j<md[ispec].nlev;j++){
           gp[gp_sorter[i]].mol[ispec].pops[j]= Ith(P,j);
         }
       }

       /* --- BROADENED FAILURE CHECK + FIXED RETRY LOGIC ---
          Original code only checked cvstatus==-3 (CV_CONV_FAILURE). Changed
          to cvstatus<0 to catch any CVode failure code, matching the more
          robust handling in the 3D solver - important since larger/stiffer
          systems (more levels) can fail in more varied ways. Also: the
          original call to reduceTol() passed reltol by value, so the
          *local* reltol in this function never actually decreased, meaning
          the surrounding while(cvstatus==-3) condition on reltol>MINTOL was
          checking a stale value and could loop indefinitely. reduceTol()
          now takes reltol by reference and updates it correctly here. */
       if(cvstatus < 0 && reltol > MINTOL){

         /* Reset P to the true initial populations before retrying,
            rather than leaving it at a possibly-diverged state. */
         for(j=0; j < NEQ; ++j) Ith(P,j) = Pops[j];

         if (using_dense == 0) {
           printf("WARNING: SPGMR failed (stiff system). Switching to Dense Direct Solver...\n");
           using_dense = 1;

           retval = CVodeInit(cvode_mem, f, radii[0], P);
           if (check_retval(&retval, "CVodeInit", 1)) printf("Failed to reinitialize CVODE\n");
           CVodeSetUserData(cvode_mem, &user_data);
           CVodeSetMaxNumSteps(cvode_mem, 5000);

           if (LS != NULL) SUNLinSolFree(LS);

           sunMatrix = SUNDenseMatrix(NEQ, NEQ);
           LS = SUNLinSol_Dense(P, sunMatrix);
           CVodeSetLinearSolver(cvode_mem, LS, sunMatrix);
           CVodeSetJacFn(cvode_mem, NULL);

         } else {
           printf("WARNING: Dense solver failed. Reducing tolerances...\n");

           retval = CVodeInit(cvode_mem, f, radii[0], P);
           if (check_retval(&retval, "CVodeInit", 1)) printf("Failed to reinitialize CVODE\n");
           CVodeSetUserData(cvode_mem, &user_data);
           CVodeSetMaxNumSteps(cvode_mem, 5000);

           if (LS != NULL) SUNLinSolFree(LS);
           if (sunMatrix != NULL) SUNMatDestroy(sunMatrix);

           sunMatrix = SUNDenseMatrix(NEQ, NEQ);
           LS = SUNLinSol_Dense(P, sunMatrix);
           CVodeSetLinearSolver(cvode_mem, LS, sunMatrix);
           CVodeSetJacFn(cvode_mem, NULL);

           reduceTol(abstol, &reltol, cvode_mem, 0.1, NEQ);
         }

         break;

       }else if(cvstatus>0){
         cvodeErrs++;
         printf("CVODE error %d (continuing to next timestep - check populations!!)\n",cvstatus);
       }

       if(cvodeErrs>=15){
         bail_out("CVODE solver failure - check physical model and tolerances.");
         exit(1);
       }
    }
  }while(cvstatus<0 && reltol > MINTOL);

  /* Free P and abstol vectors */
  N_VDestroy(P);
  N_VDestroy(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  /* Free the linear solver memory */
  if (LS != NULL) SUNLinSolFree(LS);

  /* Free the matrix memory (only allocated if we fell back to Dense) */
  if (sunMatrix != NULL) SUNMatDestroy(sunMatrix);

  free(jbar_grid);

  /* Free the heap-allocated rate matrices */
  free(p_rates);
  free(coll_rates);

  /* Free CKCdata in memory */
  if (par->useCKCdata == 1 || par->useCKCdata == 2) {
    free(CKCdata.Tedata);
    free(CKCdata.nedata);
    free(CKCdata.r_values);
  }

}


/*....................................................................*/
int
levelPops(molData *md, configInfo *par, struct grid *gp, int *popsdone, double *lamtab, double *kaptab, const int nEntries, double *radii, int *gp_sorter){

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

    freeGridCont(par, gp);
    mallocGridCont(par, md, gp);
    calcGridLinesDustOpacity(par, md, lamtab, kaptab, nEntries, gp);

    LTE(par,gp,md);

    defaultErrorHandler = gsl_set_error_handler_off();

    gridPointData *mp[par->pIntensity];
    double *halfFirstDs[par->pIntensity];

    calcGridMolSpecNumDens(par,md,gp);
    totalNMaserWarnings = 0;
    nVerticesDone=0;

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
      solveStatEq(gp,md,ispec,par,blends,nextMolWithBlend,mp,halfFirstDs, nMaserWarnings, radii, gp_sorter);
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
