/*
 *  defaults.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#include <math.h>
#include "defaults.h" /* includes lime_config.h which defines configInfo */
#include "ufunc_types.h" /* for the USERFUNC_* macros */


void
default_abundance(double x, double y, double z, double *abundance){
  abundance[0] = -1.0;
  defaultFuncFlags |= (1 << USERFUNC_abundance);
}

void
default_magfield(double x, double y, double z, double *B){
  B[0]=0.0;
  B[1]=0.0;
  B[2]=0.0;
  defaultFuncFlags |= (1 << USERFUNC_magfield);
}

void
default_gasIIdust(double x, double y, double z, double *gas2dust){
  *gas2dust=100.;
  defaultFuncFlags |= (1 << USERFUNC_gasIIdust);
}

double
default_gridDensity(configInfo *par, double *r, void (*density)(configInfo *par, double x, double y, double z, double *val)){
  /*
The grid points within the model are chosen randomly via the rejection method with a probability distribution which the present function is intended to provide.

Notes:
  - The present function is interpreted by LIME as giving *relative* probabilities, the ultimate normalization being set by the desired number of grid points conveyed to the task via par->pIntensity.
  - If par->samplingAlgorithm is chosen to be zero (the current default value), further manipulations to the probability distribution are performed according to the set value of par->sampling.
  - The user may supply their own version of the present function within model.c; the default here implements the grid-point selection function used in LIME<=1.5.
  */
  double val[99],totalDensity=0.0,rSquared=0.0,fracDensity=0.0;
  int i;

  rSquared = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
  if(rSquared>=par->radiusSqu)
    return 0.0;

  density(par,r[0],r[1],r[2],val);
  for (i=0;i<par->numDensities;i++) totalDensity += val[i];
  fracDensity = pow(totalDensity,defaultDensyPower)/par->gridDensGlobalMax;

  defaultFuncFlags |= (1 << USERFUNC_gridDensity);

  return fracDensity;
}

