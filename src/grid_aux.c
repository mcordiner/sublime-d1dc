/*
 *  grid_aux.c
 *  
 *  This file is part of SUBLIMED1D (SUBlimating gases in LIME, version D1D)
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *  Copyright (C) 2023 Martin Cordiner and Kristen Darnell (NASA GSFC)
 *
 */

#include "sublime.h"

/*....................................................................*/
void mallocAndSetDefaultGrid(struct grid **gp, const size_t numPoints, const size_t numSpecies){
  size_t i,j;

  *gp = malloc(sizeof(**gp)*numPoints);
  for(i=0;i<numPoints; i++){
    (*gp)[i].v1 = NULL;
    (*gp)[i].v2 = NULL;
    (*gp)[i].v3 = NULL;
    (*gp)[i].dir = NULL;
    (*gp)[i].neigh = NULL;
    (*gp)[i].w = NULL;
    (*gp)[i].ds = NULL;
    (*gp)[i].dens=NULL;
    (*gp)[i].t[0]=-1.0;
    (*gp)[i].t[1]=-1.0;
    (*gp)[i].B[0]=0.0;
    (*gp)[i].B[1]=0.0;
    (*gp)[i].B[2]=0.0;
    (*gp)[i].conv=0;

    if(numSpecies > 0){
      (*gp)[i].mol = malloc(sizeof(*(*gp)[i].mol)*numSpecies);
      for(j=0;j<numSpecies;j++){
        (*gp)[i].mol[j].pops        = NULL;
        (*gp)[i].mol[j].specNumDens = NULL;
        (*gp)[i].mol[j].partner     = NULL;
        (*gp)[i].mol[j].cont        = NULL;
        (*gp)[i].mol[j].dopb = 0.0;
        (*gp)[i].mol[j].binv = 0.0;
        (*gp)[i].mol[j].nmol = 0.0;
        (*gp)[i].mol[j].abun = 0.0;
      }
    }else
      (*gp)[i].mol = NULL;
  }
}

/*....................................................................*/
void calcGridMolDoppler(configInfo *par, molData *md, struct grid *gp){
  int i,id;

  for(i=0;i<par->nSpecies;i++){
    /* Calculate Doppler and thermal line broadening */
    for(id=0;id<par->ncell;id++) {
      gp[id].mol[i].dopb = sqrt(gp[id].dopb_turb*gp[id].dopb_turb\
                                + 2.*KBOLTZ/md[i].amass*gp[id].t[0]);
      gp[id].mol[i].binv = 1./gp[id].mol[i].dopb;
    }
  }
}

/*....................................................................*/
void calcGridMolDensities(configInfo *par, struct grid **gp){
  int id,ispec,i;

  for(id=0;id<par->ncell;id++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      (*gp)[id].mol[ispec].nmol = 0.0;
      for(i=0;i<par->numDensities;i++)
        (*gp)[id].mol[ispec].nmol += (*gp)[id].mol[ispec].abun*(*gp)[id].dens[i]\
                                    *par->nMolWeights[i];
    }
  }
}

/*....................................................................*/
void calcGridMolSpecNumDens(configInfo *par, molData *md, struct grid *gp){
  int gi,ispec,ei;

  for(gi=0;gi<par->ncell;gi++){
    for(ispec=0;ispec<par->nSpecies;ispec++){
      for(ei=0;ei<md[ispec].nlev;ei++){
        gp[gi].mol[ispec].specNumDens[ei] = gp[gi].mol[ispec].binv\
          *gp[gi].mol[ispec].nmol*gp[gi].mol[ispec].pops[ei];
      }
    }
  }
}


/*....................................................................*/
void distCalc(configInfo *par, struct grid *gp){
  int i,k,l;

  for(i=0;i<par->ncell;i++){
    free(gp[i].dir);
    free(gp[i].ds);
    gp[i].dir=malloc(sizeof(*(gp[i].dir)) *gp[i].numNeigh);
    gp[i].ds =malloc(sizeof(double)*gp[i].numNeigh);
    memset(gp[i].dir, 0., sizeof(*(gp[i].dir)) * gp[i].numNeigh);
    memset(gp[i].ds, 0., sizeof(double) * gp[i].numNeigh);
    for(k=0;k<gp[i].numNeigh;k++){
      for(l=0;l<3;l++)
        gp[i].dir[k].x[l] = gp[i].neigh[k]->x[l] - gp[i].x[l];

      gp[i].ds[k] = sqrt(  gp[i].dir[k].x[0]*gp[i].dir[k].x[0]\
                         + gp[i].dir[k].x[1]*gp[i].dir[k].x[1]\
                         + gp[i].dir[k].x[2]*gp[i].dir[k].x[2]);

      for(l=0;l<3;l++)
        gp[i].dir[k].xn[l] = gp[i].dir[k].x[l]/gp[i].ds[k];
    }
    gp[i].nphot=RAYS_PER_POINT;
  }
}


/*....................................................................*/
void
getEdgeVelocities(configInfo *par, struct grid *gp){
  int i,k,j,l;
  double vel[3], x[3];
  
  for(i=0;i<par->ncell;i++){
    gp[i].v1=malloc(3*gp[i].numNeigh*sizeof(double));
    gp[i].v2=malloc(3*gp[i].numNeigh*sizeof(double));
    gp[i].v3=malloc(3*gp[i].numNeigh*sizeof(double));

    for(k=0;k<gp[i].numNeigh;k++){
      for(j=0;j<3;j++) x[j]=gp[i].x[j];		
      for(l=0;l<5;l++){
        velocity(par,x[0],x[1],x[2],vel);	

        if (l==1) {
	  gp[i].v1[3*k]=vel[0]; gp[i].v1[3*k+1]=vel[1]; gp[i].v1[3*k+2]=vel[2];
        }
        if (l==2) {
          gp[i].v2[3*k]=vel[0]; gp[i].v2[3*k+1]=vel[1]; gp[i].v2[3*k+2]=vel[2];
        }
        if (l==3) {
          gp[i].v3[3*k]=vel[0]; gp[i].v3[3*k+1]=vel[1]; gp[i].v3[3*k+2]=vel[2];
        }
		
        for(j=0;j<3;j++) x[j]=x[j]+(gp[i].dir[k].xn[j]*gp[i].ds[k])/4.;
      }
    }
  }

  par->edgeVelsAvailable = 1;
}

/*....................................................................*/
int setupAndWriteGrid(configInfo *par, struct grid *gp, molData *md, char *outFileName){
  // I removed an unused function from this function
  int status = 0;
  return status;
}

/*....................................................................*/
void writeGridIfRequired(configInfo *par, struct grid *gp, molData *md, const int dataStageI){
  if(par->writeGridAtStage[dataStageI-1]){
    int status=0;
    char message[80];

    status = setupAndWriteGrid(par, gp, md, par->gridOutFiles[dataStageI-1]);

    if(status){
      sprintf(message, "writeGrid at data stage %d returned with status %d", dataStageI, status);
      if(!silent) bail_out(message);
exit(1);
    }
  }
}