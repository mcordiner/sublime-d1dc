/*
 *  lime_defaults.c
 *  SUBLIMED1D main driver script
 *  This file is part of SUBLIMED1D (SUBlimating gases in LIME, version D1D)
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *  Copyright (C) 2023 Martin Cordiner and Kristen Darnell (NASA GSFC)
 *
 */

#include "lime.h" /* includes defaults.h */

/* One of the following two must be defined by the user:  (abundance or molNumDensity)*/
void __attribute__((weak))
abundance(double x, double y, double z, double *abun){
  default_abundance(x,y,z, abun);
}

void __attribute__((weak))
magfield(double x, double y, double z, double *B){
  default_magfield(x,y,z, B);
}

void __attribute__((weak))
gasIIdust(double x, double y, double z, double *gas2dust){
  default_gasIIdust(x,y,z, gas2dust);
}

double __attribute__((weak))
gridDensity(configInfo *par, double *r){
  return default_gridDensity(par, r, density);
}

