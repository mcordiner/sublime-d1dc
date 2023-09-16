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
input(inputPars *par, image *img){
	char input[] = "input.par";
	char buffer[10000];
	char parline[10000];
	char token[1000];
	int i_str, i_dbl, i_int;
	FILE *fp;
	
	char name[200];
	char value[200];
	
	/* Strings */
	int numstr = 3;
	char *strlist[3] = {"runname", 
						"moldatfile", 
						"girdatfile"};
	char strval[numstr][100];
	
	/* doubles */
	int numdbl = 14;
	char *dbllist[14] = {"abund",
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
						"dopplerb",
						"xne",
						"tnuc"};
	double dblval[numdbl];
	
	/* ints */
	int numint = 7;
	char *intlist[7] = {"collPartIds",
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
		if (strstr(parline,"=") && strstr(parline,";")){	//Only consideres lines with the expected formatting
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
 * Basic parameters. See cheat sheet for details.
 */
 
 /* strings */
 if (strval[0] != "NaN"){
 	par->outputfile = (char*)malloc(200);
 	strcpy(par->outputfile, strval[0]);
 	strcat(par->outputfile, ".pop");    
    
    img[0].filename = (char*)malloc(200);
 	strcpy(img[0].filename, strval[0]);
 	strcat(img[0].filename, ".fits");
 }

 
 if (strval[1] != "NaN"){
 	par->moldatfile[0] = (char*)malloc(200);
 	strcpy(par->moldatfile[0], strval[1]);
 }
 
 if (strval[2] != "NaN"){
 	par->girdatfile[0] = (char*)malloc(200);
 	strcpy(par->girdatfile[0], strval[2]);
 }
  
// 	printf("par->outputfile = %s\n", par->outputfile);
// 	printf("img[0].filename = %s\n", img[0].filename);
// 	printf("par->moldatfile[0] = %s\n",par->moldatfile[0]);
// 	printf("par->girdatfile[0] = %s\n",par->girdatfile[0]);

  /* doubles */
  
  if (isnan(dblval[0]) == 0){
  	par->abund = dblval[0];
  }
  
  if (isnan(dblval[1]) == 0){
  	par->betamol = dblval[1];
  }
  
  if (isnan(dblval[2]) == 0){
  	img[0].distance = dblval[2] * AU; // source distance in m
  }
  
  if (isnan(dblval[3]) == 0){
  	img[0].imgres  = dblval[3]; // Resolution in arc seconds
  }
  
  if (isnan(dblval[4]) == 0){
  	par->Qwater  = dblval[4];
  }
  
  if (isnan(dblval[5]) == 0){
  	par->rHelio  = dblval[5]; // Heliocentric distance in AU
  }

  if (isnan(dblval[6]) == 0){
  	par->tkin  = dblval[6];
  }
  
  if (isnan(dblval[7]) == 0){
  	img[0].velres  = dblval[7]; // Channel resolution in m/s
  }
  
  if (isnan(dblval[8]) == 0){
  	par->vexp  = dblval[8];
  }
  
  if (isnan(dblval[9]) == 0){
  	par->radius  = dblval[9];
  }

  if (isnan(dblval[10]) == 0){
  	par->rnuc  = dblval[10];
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
  	par->dopplerb  = dblval[11];
  }
  
  if (isnan(dblval[12]) == 0){
  	par->xne = dblval[12];
  }

  if (isnan(dblval[13]) == 0){
  	par->tNuc  = dblval[13];
  }
  
//   	printf("par->dopplerb  = %.2e\n",par->dopplerb);
//    	printf("par->xne = %.2e\n",par->xne);
  
  /* ints */
  
  if (intval[0] != -1){
  	par->collPartIds[0]  = intval[0];
  }
  
  if (intval[1] != -1){
  	img[0].nchan  = intval[1]; // Number of channels
  }
  
  if (intval[2] != -1){
  	img[0].pxls  = intval[2]; // Pixels per dimension
  }
  
  if (intval[3] != -1){
  	img[0].trans  = intval[3]; // zero-indexed J quantum number
  }
  
  if (intval[4] != -1){
  	img[0].unit  = intval[4]; // 0:Kelvin 1:Jansky/pixel 2:SI 3:Lsun/pixel 4:tau
  }
  
//   printf("par->collPartIds[0] = %d\n",par->collPartIds[0]);
//   printf("img[0].nchan = %d\n",img[0].nchan);
//   printf("img[0].pxls = %d\n",img[0].pxls);
//   printf("img[0].trans = %d\n",img[0].trans);
//   printf("img[0].unit = %d\n",img[0].unit);
  
  /* optional ints */
  if (intval[5] != -1){
  	par->useEP  = intval[5]; // zero-indexed J quantum number
  }
  
  if (intval[6] != -1){
  	par->pIntensity = intval[6];
  }
  
//  printf("par->useEP = %d\n",par->useEP);
  
  
  
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

