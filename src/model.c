#include "lime.h"
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
	
	char name[100];
	char value[100];
	
	/* Strings */
	int numstr = 3;
	char *strlist[3] = {"runname", 
						"moldatfile", 
						"girdatfile"};
	char strval[numstr][100];
	
	/* doubles */
	int numdbl = 11;
	char *dbllist[11] = {"abund",
						"betamol",
						"delta",
						"imgres",
						"Qwater",   
						"rhelio", 
						"tkin",
						"velres",
						"vexp", 
						"dopplerb",
						"xne"};
	double dblval[numdbl];
	
	/* ints */
	int numint = 6;
	char *intlist[6] = {"collPartIds",
						"nchan", 
						"pxls", 
						"trans",
						"unit", 
						"useEP"};
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
 	par->outputfile = (char*)malloc(100);
 	strcpy(par->outputfile, "output/");
 	strcat(par->outputfile, strval[0]);
 	strcat(par->outputfile, ".pop");    
    
    img[0].filename = (char*)malloc(100);
    strcpy(img[0].filename, "output/");
 	strcat(img[0].filename, strval[0]);
 	strcat(img[0].filename, ".fits");
 }

 
 if (strval[1] != "NaN"){
 	par->moldatfile[0] = (char*)malloc(100);
 	strcpy(par->moldatfile[0], strval[1]);
 }
 
 if (strval[2] != "NaN"){
 	par->girdatfile[0] = (char*)malloc(100);
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
  	par->rHelio  = dblval[5];
  	
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
  
  if (isnan(dblval[9]) == 0){
  	par->dopplerb  = dblval[9];
  }
  
  if (isnan(dblval[10]) == 0){
  	par->xne = dblval[10];
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
  
//   printf("par->useEP = %d\n",par->useEP);
  
  
  
  /* Variables not in input.par */
  par->beta = 1.042e-5;
  par->rnuc = 2.5e2;
   
  par->radius           = 2e8;
  par->minScale         = par->rnuc;
  par->pIntensity = 500;
  par->girScale = 1.0;
  par->lte_only         = 1;
  par->useCKCdata		= 1;
  par->CKCTeFile        = "na";
  par->CKCneFile        = "na";
  par->gridfile         = "output/grid.vtk";

  par->nMolWeights[0]   = 1.0;
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

