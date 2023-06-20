#include "lime.h"
/******************************************************************************/

void
input(inputPars *par, image *img){
/*
 * Basic parameters. See cheat sheet for details.
 */
  par->beta = 1.042e-5;
  par->betamol = 1.5e-5;
  par->vexp = 700.;
  par->tkin = 50;
  par->rnuc = 2.5e2;
  par->abund = 0.001;
  par->dopplerb = 100;
   
  par->useEP       = 0;
  par->Qwater       = 1e27;
  par->xne = 0.2;
  par->rHelio       = 1.0;
  par->radius           = 2e8;
  par->minScale         = par->rnuc;
  par->pIntensity = 500;
  par->moldatfile[0]    = "data/moldat/hcn.dat";
  par->girdatfile[0]    = "data/girdat/g_hcn_1au.dat";
  par->girScale = 1.0;
  par->lte_only         = 1;
  par->useCKCdata		= 1;
  par->CKCTeFile        = "na";
  par->CKCneFile        = "na";

  par->outputfile = "output/hcn.pop";
  par->gridfile         = "output/grid.vtk";

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
  img[0].filename = "output/hcn.fits"; // Output filename
}

/******************************************************************************/

void
density(configInfo *par, double x, double y, double z, double *density){
/*
 * Define variable for radial coordinate
 */
  double r;

  const double rMin = par->rnuc; /* This cutoff should be chosen smaller than par->minScale but greater than zero (to avoid a singularity at the origin). */

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
    density[0] = par->Qwater /(4*PI*pow(r, 2)*par->vexp)*exp(-r*par->beta/par->vexp);
}

/******************************************************************************/

void
temperature(configInfo *par, double x, double y, double z, double *temperature){
  temperature[0] = par->tkin;
}

/******************************************************************************/

void
molNumDensity(configInfo *par, double x, double y, double z, double *nmol){
 /*
 * Define variable for radial coordinate
 */
  double r;

  const double rMin = par->rnuc; /* This cutoff should be chosen smaller than par->minScale but greater than zero (to avoid a singularity at the origin). */

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
    nmol[0] =par->abund*par->Qwater/(4*PI*pow(r, 2)*par->vexp)*exp(-r*par->betamol/par->vexp);
}

/******************************************************************************/

void
doppler(configInfo *par, double x, double y, double z, double *doppler){
  *doppler = par->dopplerb;
}

/******************************************************************************/

void
velocity(configInfo *par, double x, double y, double z, double *vel){
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
  vel[0]=par->vexp*sin(theta)*cos(phi);
  vel[1]=par->vexp*sin(theta)*sin(phi);
  vel[2]=par->vexp*cos(theta);
}

