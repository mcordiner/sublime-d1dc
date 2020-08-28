#include "lime.h"

double beta= 1.042e-5;
double betahcn= 1.5e-5;

double vexp = 700.;
double Qwater = 1e27;
double tkin = 50.;
double rnuc = 2.5e2;
double abund = 0.001;

/******************************************************************************/

void
input(inputPars *par, image *img){
/*
 * Basic parameters. See cheat sheet for details.
 */
  par-> useEP       = 0;
  par->Qwater       = Qwater;
  par -> xne = 0.2;
  par->rHelio       = 1.0;
  par->radius           = 1e8;
  par->minScale         = rnuc;
  par->pIntensity = 1000;
  par->sinkPoints = 500;
  par->moldatfile[0]    = "hcn.dat";
  par->girdatfile[0]    = "g_hcn_1au.dat";
  par->lte_only         = 0;

  par->nSolveIters  = 7;
  par->outputfile = "hcn.pop";
  par->binoutputfile    = "restart.pop";
  par->gridfile         = "grid.vtk";

  par->collPartIds[0]   = 1;
  par->nMolWeights[0]   = 1.0;

/*
 * Definitions for image #0. Add blocks for additional images.
 */
  img[0].velres         = 100.;   // Channel resolution in m/s
  img[0].nchan          = 50;     // Number of channels
  img[0].trans          = 3;    // zero-indexed J quantum number
  img[0].pxls           = 256;    // Pixels per dimension
  img[0].imgres         = 0.5;    // Resolution in arc seconds
  img[0].distance = 1.0*AU; // source distance in m
  img[0].unit           = 0;    // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  img[0].filename = "hcn.fits"; // Output filename
}

/******************************************************************************/

void
density(double x, double y, double z, double *density){
/*
 * Define variable for radial coordinate
 */
  double r;

  const double rMin = rnuc; /* This cutoff should be chosen smaller than par->minScale but greater than zero (to avoid a singularity at the origin). */

  /*
   * Calculate radial distance from origin
   */
  r=sqrt(x*x+y*y+z*z);
  /*
   * Calculate a Haser density profile
   * (Multiply with 1e6 to go to SI-units)
   */
  if(r<rMin)
    density[0] = 1e-20; /* Just to prevent overflows at r==0! */
  else
    density[0] = Qwater /(4*PI*pow(r, 2)*vexp)*exp(-r*beta/vexp);
}

/******************************************************************************/

void
temperature(double x, double y, double z, double *temperature){
  temperature[0] = tkin;
}

/******************************************************************************/

void
molNumDensity(double x, double y, double z, double *nmol){
 /*
 * Define variable for radial coordinate
 */
  double r;

  const double rMin = rnuc; /* This cutoff should be chosen smaller than par->minScale but greater than zero (to avoid a singularity at the origin). */

  /*
   * Calculate radial distance from origin
   */
  r=sqrt(x*x+y*y+z*z);
  /*
   * Calculate a Haser density profile
   * (Multiply with 1e6 to go to SI-units)
   */
  if(r<rMin)
    nmol[0] = 0.;
  else
    nmol[0] =abund*Qwater/(4*PI*pow(r, 2)*vexp)*exp(-r*betahcn/vexp);
}

/******************************************************************************/

void
doppler(double x, double y, double z, double *doppler){
  *doppler = 100.;
}

/******************************************************************************/

void
velocity(double x, double y, double z, double *vel){
/*
 * Variables for spherical coordinates
 */
  double phi, theta;
/*
 * Transform Cartesian coordinates into spherical coordinates
 */
  theta=atan2(sqrt(x*x+y*y),z);
  phi=atan2(y,x);
/*
 * Vector transformation back into Cartesian basis
 */
  vel[0]=vexp*sin(theta)*cos(phi);
  vel[1]=vexp*sin(theta)*sin(phi);
  vel[2]=vexp*cos(theta);
}
