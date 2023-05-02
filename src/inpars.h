/*
 *  inpars.h
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#ifndef INPARS_H
#define INPARS_H

#include "dims.h"

/* input parameters */
typedef struct {
  double beta,betamol,vexp,tkin,rnuc,abund; // Formerly magic numbers in model.c
  double radius,minScale,tcmb,Qwater,rHelio,xne,colliScale,girScale,tNuc,*nMolWeights,*dustWeights;
  double (*gridDensMaxLoc)[DIM],*gridDensMaxValues,*collPartMolWeights;
  int sinkPoints,pIntensity,blend,*collPartIds,traceRayAlgorithm,samplingAlgorithm;
  int sampling,lte_only,init_lte,antialias,polarization,nThreads,nSolveIters,useEP,fixRNG,useCKCdata;
  char **girdatfile,**moldatfile,**collPartNames;
  char *outputfile,*binoutputfile,*gridfile,*pregrid,*restart,*dust,*CKCTeFile,*CKCneFile;
  char *gridInFile,**gridOutFiles;
  _Bool resetRNG,doSolveRTE;
} inputPars;

/* Image information */
typedef struct {
  int nchan,trans,molI;
  double velres;
  double imgres;
  int pxls;
  int unit;
  char *units;
  int psfShape;
  double psfWidth;
  int rebinSpec;
  int nBins;
  double binWidth;
  double freq,bandwidth;
  char *filename;
  double source_vel;
  double theta,phi,incl,posang,azimuth;
  double distance;
  _Bool doInterpolateVels;
} image;

#endif /* INPARS_H */
