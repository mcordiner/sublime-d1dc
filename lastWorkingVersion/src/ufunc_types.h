/*
 *  ufunc_types.h
 *  SUBLIMED1D main driver script
 *  This file is part of SUBLIMED1D (SUBlimating gases in LIME, version D1D)
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *  Copyright (C) 2023 Martin Cordiner and Kristen Darnell (NASA GSFC)
 *
 */


#ifndef UFUNC_TYPES_H
#define UFUNC_TYPES_H

#include "lime_config.h" /* for configInfo */

#define USERFUNC_density       0
#define USERFUNC_temperature   1
#define USERFUNC_abundance     2
#define USERFUNC_molNumDensity 3
#define USERFUNC_doppler       4
#define USERFUNC_velocity      5
#define USERFUNC_magfield      6
#define USERFUNC_gasIIdust     7
#define USERFUNC_gridDensity   8

//extern int defaultFuncFlags;

void	density(configInfo*,double,double,double,double *);
void	temperature(configInfo*,double,double,double,double *);
void	abundance(double,double,double,double *);
void	molNumDensity(configInfo*,double,double,double,double *);
void	doppler(configInfo*,double,double,double, double *);
void	velocity(configInfo*,double,double,double,double *);
void	magfield(double,double,double,double *);
void	gasIIdust(double,double,double,double *);
double	gridDensity(configInfo*, double*);

#endif /* UFUNC_TYPES_H */

