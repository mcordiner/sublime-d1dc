/*
 *  defaults.h
 *  
 *  This file is part of SUBLIMED1D (SUBlimating gases in LIME, version D1D)
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *  Copyright (C) 2023 Martin Cordiner and Kristen Darnell (NASA GSFC)
 *

Files which include this must
  - define a struct named configInfo which has members as follows:
      double radiusSqu;
      double gridDensGlobalMax;
      int numDensities;

  - define macros USERFUNC_* as used in defaults.c.
 */

#ifndef DEFAULTS_H
#define DEFAULTS_H

#include "lime_config.h" /* for configInfo. */

#define DENSITY_POWER           0.2

extern int defaultFuncFlags;
extern double defaultDensyPower;

void	default_abundance(double x, double y, double z, double *abundance);
void	default_magfield(double x, double y, double z, double *B);
void	default_gasIIdust(double x, double y, double z, double *gas2dust);
double	default_gridDensity(configInfo *par, double *r, void (*density)(configInfo *par, double x, double y, double z, double *val));

#endif /* DEFAULTS_H */

