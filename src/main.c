/*
 *  main.c
 *  
 *  This file is part of SUBLIMED1D (SUBlimating gases in LIME, version D1D)
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *  Copyright (C) 2023 Martin Cordiner and Kristen Darnell (NASA GSFC)
 *
 */


#include "sublime.h"

#ifdef NOVERBOSE
int silent = 1;
#else
int silent = 0;
#endif

/*....................................................................*/
int
initParImg(inputPars *par, image **img)
{
  /* Initialize par with default values, allocate space for the
     output fits images, initialize the images with default values,
     and finally call the input() routine from model.c to set both the
     par and image values.
  */

  int i,j,id,nImages;
  const double defaultAngle=-999.0;

  /* Set 'impossible' default values for mandatory parameters */
  par->radius    = 0;
  par->minScale  = 0;
  par->sinkPoints= 0;
  par->Qwater    = 0;
  par->rHelio    = 0; 
  par->beta 	 = 0;
  par->vexp 	 = 0;
  par->tkin 	 = 0;
  par->rnuc 	 = 0;
  par->abund 	 = 0;

  /* Set default values for optional parameters */
  par->dust  	    = NULL;
  par->outputfile   = NULL;
  par->binoutputfile= NULL;
  par->gridfile     = NULL;
  par->pregrid      = NULL;
  par->restart      = NULL;
  par->gridInFile   = NULL;
  par->CKCTeFile    = NULL;
  par->CKCneFile    = NULL;
  
  par->collPartIds  = malloc(sizeof(int)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->collPartIds[i] = 0; /* Possible values start at 1. */
  par->nMolWeights  = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->nMolWeights[i] = -1.0;
  par->dustWeights  = malloc(sizeof(double)*MAX_N_COLL_PART); /* This param no longer has any effect. */
  for(i=0;i<MAX_N_COLL_PART;i++) par->dustWeights[i] = -1.0;
  par->collPartMolWeights = malloc(sizeof(double)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++) par->collPartMolWeights[i] = -1.0;

  par->gridDensMaxValues = malloc(sizeof(*(par->gridDensMaxValues))*MAX_N_HIGH);
  par->gridDensMaxLoc    = malloc(sizeof(*(par->gridDensMaxLoc))*MAX_N_HIGH);
  for(i=0;i<MAX_N_HIGH;i++){
    par->gridDensMaxValues[i] = -1.0; /* Impossible default value. */
    for(j=0;j<DIM;j++) par->gridDensMaxLoc[i][j] = 0.0;
  }

  par->tcmb = LOCAL_CMB_TEMP;
  par->lte_only=0;
  par->init_lte=0;
  par->samplingAlgorithm=0;
  par->sampling=2; /* Now only accessed if par->samplingAlgorithm==0. */
  par->blend=0;
  par->antialias=1;
  par->polarization=0;
  par->nThreads = NTHREADS;
  par->nSolveIters=1;
  par->traceRayAlgorithm=0;
  par->resetRNG=0;
  par->doSolveRTE=0;
  par->colliScale = 1.0;
  par->girScale = 1.0;
  par->xne = DEFAULT_XNE;
  par->useEP = 0;
  par->fixRNG = 1;
  par->tNuc = DEFAULT_TNUC;
  par->dopplerb=0;
  par->useCKCdata = 0;
  par->beta = DEFAULT_BETA;
  par->pIntensity= DEFAULT_PINTENSITY;

  par->gridOutFiles = malloc(sizeof(char *)*NUM_GRID_STAGES);
  for(i=0;i<NUM_GRID_STAGES;i++)
    par->gridOutFiles[i] = NULL;

  /* Allocate initial space for molecular data file names */
  par->moldatfile=malloc(sizeof(char *)*MAX_NSPECIES);
  par->girdatfile=malloc(sizeof(char *)*MAX_NSPECIES);
  for(id=0;id<MAX_NSPECIES;id++){
    par->moldatfile[id]=NULL;
    par->girdatfile[id]=NULL;
  }

  /* Allocate initial space for (non-LAMDA) collision partner names */
  par->collPartNames=malloc(sizeof(char *)*MAX_N_COLL_PART);
  for(i=0;i<MAX_N_COLL_PART;i++){
    par->collPartNames[i]=NULL;
  }

  /* Allocate initial space for output fits images */
  (*img)=malloc(sizeof(**img)*MAX_NIMAGES);
  for(i=0;i<MAX_NIMAGES;i++){
    (*img)[i].filename=NULL;
    (*img)[i].units=NULL;
  }

  /* First call to the user function which sets par, img values. Note that, as far as img is concerned, here we just want to find out how many images the user wants, so we can malloc the array properly. After setting up the default image parameters below, we call input() a second time then to set the actual per-image parameter values. SUBLIMED1D only makes one image so two calls are not required. If multiple images are implemented in the future, then we need to change how the beta heliocentric scaling is done inside input() or else the beta values will be scaled twice!
  */
  //input(par, *img);

  /* If the user has provided a list of image filenames, the corresponding elements of (*img).filename will be non-NULL. Thus we can deduce the number of images from the number of non-NULL elements. */
  //nImages=0;
  //while((*img)[nImages].filename!=NULL && nImages<MAX_NIMAGES)
  // nImages++;

   nImages=1;

  /* Set img defaults. */
  for(i=0;i<nImages;i++) {
    (*img)[i].nchan      =  0;
    (*img)[i].trans      = -1;
    (*img)[i].molI       = -1;
    (*img)[i].velres     = -1.0;
    (*img)[i].imgres     = -1.0;
    (*img)[i].pxls       = -1;
    (*img)[i].unit       =  0;
    (*img)[i].freq       = -1.0;
    (*img)[i].bandwidth  = -1.0;
    (*img)[i].psfShape   = 0;
    (*img)[i].psfWidth   = 0.0;
    (*img)[i].rebinSpec   = 0;
    (*img)[i].nBins      = 0;
    (*img)[i].binWidth   = 0;
    (*img)[i].source_vel =  0.0;
    (*img)[i].theta      =  0.0;
    (*img)[i].phi        =  0.0;
    (*img)[i].incl       = defaultAngle;
    (*img)[i].posang     = defaultAngle;
    (*img)[i].azimuth    = defaultAngle;
    (*img)[i].distance   = -1.0;
    (*img)[i].doInterpolateVels = FALSE;
  }

  /* LIME requests a second-pass reading of the user-set parameters (this time just to read the par->moldatfile and img stuff), but presently, we are only making a single image so don't need to. */
  input(par,*img);

  return nImages;


}

/*....................................................................*/
int main() {
  /* Main program for stand-alone SUBLIME */

  inputPars par;
  image *img = NULL;
  int nImages,status=0;
  char message[STR_LEN_0];

  (void)status; // just to stop compiler warnings because this return value is currently unused.

  nImages = initParImg(&par, &img);

  status = run(par, img, nImages);
  if(status){
    if(!silent){
      sprintf(message, "Function run() returned with status %d", status);
      bail_out(message);
    }
exit(1);
  }

  free(img);
  freeInputPars(&par);

  return 0;
}

