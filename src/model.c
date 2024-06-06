/*
 *  model.c
 *
 *  This file is part of SUBLIMED1D (SUBlimating gases in LIME, version D1D)
 *
 *  Copyright (C) 2023 Martin Cordiner and Kristen Darnell (NASA GSFC)
 *
 */


#include "sublime.h"
#include <ctype.h>
/******************************************************************************/

void 
remove_spaces (char* restrict str_trimmed, const char* restrict str_untrimmed) // From: https://stackoverflow.com/questions/1726302/remove-spaces-from-a-string-in-c
{
  while (*str_untrimmed != '\0')
  {
    if(!isspace(*str_untrimmed))
    {
      *str_trimmed = *str_untrimmed;
      str_trimmed++;
    }
    str_untrimmed++;
  }
  *str_trimmed = '\0';
}

int 
indexOf(char *elm, char *arr[], int arr_cnt)
{
    // decreasing array count till it reaches negative
    // arr_cnt - 1 to 0
    while (arr_cnt--)
    {
        // Return array index if current element equals provided element
        if (strcmp(arr[arr_cnt],elm)==0)
            return arr_cnt;
    }

    // Element not present
    return -1; // Should never reaches this point
}

void
input(inputPars *par, image *img, char *input){
	char buffer[1000];
	char parline[1000];
	char token[1000];
	int i_str, i_dbl, i_int;
	FILE *fp;
	
	char name[200];
	char value[200];
	
	/* Strings */
	int numstr = 5;
	char *strlist[5] = {"runname", 
						"moldatfile", 
						"girdatfile",
						"moldatfile2",
						"girdatfile2"
						};
	char strval[numstr][1000];
	
	/* doubles */
	int numdbl = 22;
	char *dbllist[22] = {"abund",
						"betamol",
						"delta",
						"imgres",
						"Qwater",   
						"rhelio", 
						"tkin",
						"velres",
						"vexp",
						"radius",
						"rnuc",
						"lp",
						"dAbund", 
						"dopplerb",
						"xne",
						"tnuc",
						"colliScale",
						"ratio",
						"freq",
						"deldot",
						"xsec",
						"fwhm"};
	double dblval[numdbl];
	
	/* ints */
	int numint = 7;
	char *intlist[7] = {"collPartId",
						"nchan", 
						"pxls", 
						"trans",
						"unit", 
						"useEP",
						"npts"};
	int intval[numint];

	/* Populate arrays */
	for(int i_str = 0; i_str < numstr; i_str++) {
		strcpy(strval[i_str],"NaN");
	}
	for(int i_dbl = 0; i_dbl < numdbl; i_dbl++) {
		dblval[i_dbl] = NAN;
	} 
	for(int i_int = 0; i_int < numint; i_int++) {
		intval[i_int] = -1;
	} 
	
	/* Read input file and assign data to arrays */
	fp = fopen (input, "r");

	while(fgets(buffer, sizeof buffer, fp) != NULL) {
		remove_spaces(parline, buffer);	// Removes white space from string
		if (strstr(parline,"=") && strstr(parline,";")){	//Only considers lines with the expected formatting
			strcpy(name,strtok(parline, "="));
			strcpy(value,strtok(NULL, ";"));			
			i_str = indexOf(name,strlist,numstr);
			i_dbl = indexOf(name,dbllist,numdbl);
			i_int = indexOf(name,intlist,numint);
			if (i_str != -1){
				strcpy(strval[i_str],value);
			}
			if (i_dbl != -1){
				dblval[i_dbl] = atof(value);
			}
			if (i_int != -1){
				intval[i_int] = atoi(value);
			}
		}
		
	}
	fclose(fp);
	
/*
 * Basic parameters.
 */
 
 /* strings */
 if (strcmp(strval[0],"NaN")){
 	par->outputfile = (char*)malloc(1000);
 	strcpy(par->outputfile, strval[0]);
 	strcat(par->outputfile, ".pop");    
    
    img[0].filename = (char*)malloc(1000);
 	strcpy(img[0].filename, strval[0]);
 	strcat(img[0].filename, ".fits");
 }else{
   bail_out("Required parameter 'runname' not read correctly from input.par");
   exit(1);
  }

 
 if (strcmp(strval[1],"NaN")){
 	par->moldatfile[0] = (char*)malloc(1000);
 	strcpy(par->moldatfile[0], strval[1]);
 }else{
   bail_out("Required parameter 'moldatfile' not read correctly from input.par");
   exit(1);
  }
 
 if (strcmp(strval[2],"NaN")){
 	par->girdatfile[0] = (char*)malloc(1000);
 	strcpy(par->girdatfile[0], strval[2]);
 }
  
// 	printf("par->outputfile = %s\n", par->outputfile);
// 	printf("img[0].filename = %s\n", img[0].filename);
// 	printf("par->moldatfile[0] = %s\n",par->moldatfile[0]);
// 	printf("par->girdatfile[0] = %s\n",par->girdatfile[0]);

  /* doubles */
  
  if (isnan(dblval[0]) == 0){
  	par->abund = dblval[0];
  }else{
   bail_out("Required parameter 'abund' not read correctly from input.par");
   exit(1);
  }
  
  if (isnan(dblval[1]) == 0){
  	par->betamol = dblval[1];
  }else{
   bail_out("Required parameter 'betamol' not read correctly from input.par");
   exit(1);
  }
  
  if (isnan(dblval[2]) == 0){
  	img[0].distance = dblval[2] * AU; // source distance in m
  }else{
   bail_out("Required parameter 'delta' not read correctly from input.par");
   exit(1);
  }
  
  if (isnan(dblval[3]) == 0){
  	img[0].imgres  = dblval[3]; // Resolution in arc seconds
  }else{
   bail_out("Required parameter 'imgres' not read correctly from input.par");
   exit(1);
  }
  
  if (isnan(dblval[4]) == 0){
  	par->Qwater  = dblval[4];
  }else{
   bail_out("Required parameter 'Qwater' not read correctly from input.par");
   exit(1);
  }
  
  if (isnan(dblval[5]) == 0){
  	par->rHelio  = dblval[5]; // Heliocentric distance in AU
  }else{
   bail_out("Required parameter 'rhelio' not read correctly from input.par");
   exit(1);
  }

  if (isnan(dblval[6]) == 0){
  	par->tkin  = dblval[6];
  }else{
   bail_out("Required parameter 'tkin' not read correctly from input.par");
   exit(1);
  }
  
  if (isnan(dblval[7]) == 0){
  	img[0].velres  = dblval[7]; // Channel resolution in m/s
  }else{
   bail_out("Required parameter 'velres' not read correctly from input.par");
   exit(1);
  }
  
  if (isnan(dblval[8]) == 0){
  	par->vexp  = dblval[8];
  }else{
   bail_out("Required parameter 'vexp' not read correctly from input.par");
   exit(1);
  }
  
  if (isnan(dblval[9]) == 0){
  	par->radius  = dblval[9];
  }else{
   bail_out("Required parameter 'radius' not read correctly from input.par");
   exit(1);
  }

  if (isnan(dblval[10]) == 0){
  	par->rnuc  = dblval[10];
  }else{
   bail_out("Required parameter 'rnuc' not read correctly from input.par");
   exit(1);
  }

 
//   printf("par->abund  = %.2e\n",par->abund);
//   printf("par->betamol  = %.2e\n",par->betamol);
//   printf("img[0].distance  = %.2e\n",img[0].distance);
//   printf("img[0].imgres = %.2e\n",img[0].imgres);
//   printf("par->Qwater  = %.2e\n",par->Qwater);
//   printf("par->rHelio  = %.2e\n",par->rHelio);
//   printf("par->tkin  = %.2e\n",par->tkin);
//   printf("img[0].velres  = %.2e\n",img[0].velres);
//   printf("par->vexp  = %.2e\n",par->vexp);
  
  
  /* optional doubles */

  if (isnan(dblval[11]) == 0){
  	par->lp  = dblval[11];
  }
  
  if (isnan(dblval[12]) == 0){
  	par->dAbund  = dblval[12];
  }
  
  if (isnan(dblval[13]) == 0){
  	par->dopplerb  = dblval[13];
  }
  
  if (isnan(dblval[14]) == 0){
  	par->xne = dblval[14];
  }

  if (isnan(dblval[15]) == 0){
  	par->tNuc  = dblval[15];
  }
  
  if (isnan(dblval[16]) == 0){
  	par->colliScale  = dblval[16];
  }
  
  if (isnan(dblval[17]) == 0){
  	par->ratio  = dblval[17];
  }
  
  if (isnan(dblval[18]) == 0){
  	img[0].freq  = dblval[18];
  }
  
  if (isnan(dblval[19]) == 0){
  	img[0].source_vel  = dblval[19];
  }
  
  if (isnan(dblval[20]) == 0){
  	par->xsec  = dblval[20];
  }
  
  if (isnan(dblval[21]) == 0){
  	img[0].fwhm  = dblval[21];
  }
  
  /* ints */
  
  if (intval[0] != -1){
  	par->collPartIds[0]  = intval[0];
  }else{
   bail_out("Required parameter 'collPartId' not read correctly from input.par");
   exit(1);
  }
  
  if (intval[1] != -1){
  	img[0].nchan  = intval[1]; // Number of channels
  }else{
   bail_out("Required parameter 'nchan' not read correctly from input.par");
   exit(1);
  }
  
  if (intval[2] != -1){
  	img[0].pxls  = intval[2]; // Pixels per dimension
  }else{
   bail_out("Required parameter 'pxls' not read correctly from input.par");
   exit(1);
  }
  
  if (intval[4] != -1){
  	img[0].unit  = intval[4]; // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  }else{
   bail_out("Required parameter 'unit' not read correctly from input.par");
   exit(1);
  }
  
//   printf("par->collPartIds[0] = %d\n",par->collPartIds[0]);
//   printf("img[0].nchan = %d\n",img[0].nchan);
//   printf("img[0].pxls = %d\n",img[0].pxls);
//   printf("img[0].trans = %d\n",img[0].trans);
//   printf("img[0].unit = %d\n",img[0].unit);
  
  /* optional ints */
  if (intval[3] != -1){
  	img[0].trans  = intval[3]; // zero-indexed J quantum number
  }
  
  if (intval[5] != -1){
  	par->useEP  = intval[5]; // zero-indexed J quantum number
  }
  
  if (intval[6] != -1){
  	par->pIntensity = intval[6];
  }
  
//  printf("par->useEP = %d\n",par->useEP);
  
  /* optional strings */
  if (strcmp(strval[3],"NaN")){
 	par->moldatfile[1] = (char*)malloc(1000);
 	strcpy(par->moldatfile[1], strval[3]);
  }
  if (strcmp(strval[4],"NaN")){
 	par->girdatfile[1] = (char*)malloc(1000);
 	strcpy(par->girdatfile[1], strval[4]);
  }
  
  
  /* Other model variables the user might want to change */
  
  par->minScale         = par->rnuc;
  par->girScale = 1.0/pow(par->rHelio,2);
  par->lte_only         = 0;
  par->useCKCdata		= 0;
  par->CKCTeFile        = "na";
  par->CKCneFile        = "na";
  par->gridfile         = "grid.vtk";

  par->nMolWeights[0]   = 1.0;

  /* Rescale photo-rates by heliocentric distance */
  par->beta = par->beta * par->girScale;
  par->betamol = par->betamol * par->girScale;
  
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
   * Calculate Haser density profile
   */
  if(r<rMin)
    density[0] = 1e-20; /* Prevent overflows at r==0 */
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
  double r,ld;

  const double rMin = par->rnuc; /* This cutoff should be chosen smaller than par->minScale but greater than zero (to avoid a singularity at the origin). */
  
  /*
   * Haser daughter scale length
   */
  ld = par->vexp/par->betamol;

  /*
   * Calculate radial distance from origin
   */
  r=sqrt(x*x+y*y+z*z);
  /*
   * Calculate Haser parent + daughter density profile
   * Parent and daughter components can be independently specified (or set to zero).  
   * Two molecular species can be specified with an abundance ratio given by par->ratio (e.g. A/E CH3OH).
   * If the lp value is 0, we end up with the strangely redundant case of a dual parent distribution,
   * so best to set input dAbund to zero in that case, to avoid unwanted behaviour. 
   */
  
  if(par->lp > 0.){ // Avoid dividing by zero in case of zero parent scale length
  
  if(par->ratio == 0.0){
    if(r<rMin){
      nmol[0] = 0.;
    }else{
      nmol[0] =par->abund*par->Qwater/(4*PI*pow(r, 2)*par->vexp)*exp(-r*par->betamol/par->vexp) + (par->dAbund * par->Qwater * (ld/(par->lp-ld)) * (exp(-r / par->lp) - exp(-r / ld)) / (4.0 * PI * pow(r, 2)*par->vexp));
    }
  }else{  
    if(r<rMin){
      nmol[0] = 0.;
      nmol[1] = 0.;
    }else{
      nmol[0] = par->abund * (par->ratio/(1+par->ratio)) * par->Qwater/(4*PI*pow(r, 2)*par->vexp)*exp(-r*par->betamol/par->vexp) + (par->dAbund * (par->ratio/(1+par->ratio)) * par->Qwater * (ld/(par->lp-ld)) * (exp(-r / par->lp) - exp(-r / ld)) / (4.0 * PI * pow(r, 2)*par->vexp));
      nmol[1] = par->abund * (1/(1+par->ratio)) * par->Qwater/(4*PI*pow(r, 2)*par->vexp)*exp(-r*par->betamol/par->vexp) + (par->dAbund * (1/(1+par->ratio)) * par->Qwater * (ld/(par->lp-ld)) * (exp(-r / par->lp) - exp(-r / ld)) / (4.0 * PI * pow(r, 2)*par->vexp));
    }
  }
  
  }else{
  
   if(par->ratio == 0.0){
    if(r<rMin){
      nmol[0] = 0.;
    }else{
      nmol[0] = par->abund*par->Qwater/(4*PI*pow(r, 2)*par->vexp)*exp(-r*par->betamol/par->vexp) + (par->dAbund * par->Qwater/(4*PI*pow(r, 2)*par->vexp)*exp(-r*par->betamol/par->vexp));
    }
  }else{  
    if(r<rMin){
      nmol[0] = 0.;
      nmol[1] = 0.;
    }else{
      nmol[0] = par->abund * (par->ratio/(1+par->ratio)) * par->Qwater/(4*PI*pow(r, 2)*par->vexp)*exp(-r*par->betamol/par->vexp) + (par->dAbund * (par->ratio/(1+par->ratio)) * par->Qwater/(4*PI*pow(r, 2)*par->vexp)*exp(-r*par->betamol/par->vexp));
      nmol[1] = par->abund * (1/(1+par->ratio)) * par->Qwater/(4*PI*pow(r, 2)*par->vexp)*exp(-r*par->betamol/par->vexp) + (par->dAbund * (1/(1+par->ratio)) * par->Qwater/(4*PI*pow(r, 2)*par->vexp)*exp(-r*par->betamol/par->vexp));
    }
  } 
  
  }
  
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

