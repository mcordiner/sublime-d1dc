/*
 *  raytrace.c
 *  This file is part of LIME, the versatile line modeling engine
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *
 */

#include "lime.h"
#include "raythrucells.h"
#include <math.h>


struct flux{
  double *intense;
  double *tau;
};

struct radius_struct{
  int *id;
  double *radius;
  struct flux *flux;
};

// double linear_interp(double x_array, double y_array,configInfo *par, double value){

//   gsl_spline *spline = NULL;
//   gsl_interp_accel *acc = NULL;
//   double *knus=NULL, *dusts=NULL;

//   acc = gsl_interp_accel_alloc();
//   spline = gsl_spline_alloc(gsl_interp_linear,par->pIntensity);
//   gsl_spline_init(spline,x_array,y_array,par->pIntensity);

//   interp_value = gsl_spline_eval (spline, value, acc);

//   gsl_spline_free(spline);
//   gsl_interp_accel_free(acc);

//   return(interp_value);
// }

double linear_interp(double x0, double x1, double y0, double y1, double value){
  return((y0*(x1-value) + y1*(value-x0))/(x1-x0));
}

/*....................................................................*/
void calcGridContDustOpacity(configInfo *par, const double freq\
  , double *lamtab, double *kaptab, const int nEntries, struct grid *gp){

  int id;
  double gtd;
  gsl_spline *spline = NULL;
  gsl_interp_accel *acc = NULL;
  double *kappatab = NULL;
  double *knus=NULL, *dusts=NULL;
  double *freqs=NULL;

  kappatab = malloc(sizeof(*kappatab)*1);
  knus     = malloc(sizeof(*knus)    *1);
  dusts    = malloc(sizeof(*dusts)   *1);
  freqs    = malloc(sizeof(*freqs)   *1);

  freqs[0] = freq;

  if(par->dust == NULL)
    kappatab[0] = 0.;
  else{
    acc = gsl_interp_accel_alloc();
    spline = gsl_spline_alloc(gsl_interp_cspline,nEntries);
    gsl_spline_init(spline,lamtab,kaptab,nEntries);
    kappatab[0] = interpolateKappa(freq, lamtab, kaptab, nEntries, spline, acc);
  }

  for(id=0;id<par->ncell;id++){
    gasIIdust(gp[id].x[0],gp[id].x[1],gp[id].x[2],&gtd);
    calcDustData(par, gp[id].dens, freqs, gtd, kappatab, 1, gp[id].t, knus, dusts); /* in aux.c. */
    gp[id].cont.knu = knus[0];
    gp[id].cont.dust = dusts[0];
  }

  if(par->dust != NULL){
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  free(knus);
  free(dusts);
  free(freqs);
  free(kappatab);
}


/*....................................................................*/
void
traceray(imageInfo *img,configInfo *par,struct grid *gp,molData *md,struct radius_struct radius_struct,const int im, int index){

  int ichan,stokesId,di,i,posn,posnu,posnl,molI,lineI;
  double zp,x[DIM],dx[DIM],dz,dtau,col,r;
  double contJnu,contAlpha,jnu,alpha,lineRedShift,vThisChan,deltav,vfac=0.0;
  double remnantSnu,expDTau,brightnessIncrement;

  zp=-par->radiusSqu;

  for(di=0;di<DIM;di++){
    dx[di]= img[im].rotMat[di][2]; /* This points away from the observer. */
  }
  x[0] = radius_struct.radius[index];
  x[1] = 0.0;
  x[2] = zp;

  col=0.0;
  do{
    r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    if (fabs(r) > par->minScale)
      dz = r / 5.0;
    else
      dz = 2*par->minScale;

    if (r > radius_struct.radius[0] && r < radius_struct.radius[par->pIntensity - 1]){
      for (i = 0; i < par->pIntensity; i++){
        if (radius_struct.radius[i] > r) {
          if(i==0){
            posn = radius_struct.id[i];
            dz = fabs(radius_struct.radius[i]-r);
            break;
          }
          else{
            if(fabs(radius_struct.radius[i]-r) < fabs(r-radius_struct.radius[i-1])){
              posn = radius_struct.id[i];
              break;
            }

            else{
              posn = radius_struct.id[i-1];
              break;
            }
          }//end else
        }// end if (radius_struct.radius[i] > r)
      }

      /* Calculate first the continuum stuff because it is the same for all channels:*/
      contJnu = 0.0;
      contAlpha = 0.0; 
      sourceFunc_cont(gp[posn].cont, &contJnu, &contAlpha);

      for(ichan=0;ichan<img[im].nchan;ichan++){
        jnu = contJnu;
        alpha = contAlpha;
        vThisChan = (ichan-(img[im].nchan-1)*0.5)*img[im].velres; /* Consistent with the WCS definition in writefits(). */

        if(img[im].doline){
          for(molI=0;molI<par->nSpecies;molI++){
            for(lineI=0;lineI<md[molI].nline;lineI++){
              if(md[molI].freq[lineI] > img[im].freq-img[im].bandwidth*0.5\
              && md[molI].freq[lineI] < img[im].freq+img[im].bandwidth*0.5){
                /* Calculate the red shift of the transition wrt to the frequency specified for the image.
                */
                if(img[im].trans > -1){
                  lineRedShift=(md[molI].freq[img[im].trans]-md[molI].freq[lineI])/md[molI].freq[img[im].trans]*CLIGHT;
                } else {
                  lineRedShift=(img[im].freq-md[molI].freq[lineI])/img[im].freq*CLIGHT;
                }
                deltav = vThisChan - img[im].source_vel - lineRedShift;

                velocity(x[0],x[1],x[2],gp[posn].vel);
                vfac = gaussline(deltav-dotProduct3D(dx,gp[posn].vel),gp[posn].mol[molI].binv);

                /* Increment jnu and alpha for this Voronoi cell by the amounts appropriate to the spectral line. */
                sourceFunc_line(&md[molI],vfac,&(gp[posn].mol[molI]),lineI,&jnu,&alpha);

              }//end if freq and bandwidth
            }//end for lineI
          }//end for molI
        }//end if doline
        dtau=alpha*dz;
        calcSourceFn(dtau, par, &remnantSnu, &expDTau);
        remnantSnu *= jnu*dz;

#ifdef FASTEXP
        brightnessIncrement = FastExp(radius_struct.flux[index].tau[ichan])*remnantSnu;
#else
        brightnessIncrement =    exp(-radius_struct.flux[index].tau[ichan])*remnantSnu;
#endif

        radius_struct.flux[index].intense[ichan] += brightnessIncrement;
        radius_struct.flux[index].tau[ichan] += dtau;

      }//end for ichan
    }// end if
    col+=dz;
    x[2] += dz;
    
  }while(col < 2.0*fabs(zp));
}

/*....................................................................*/
void
raytrace(int im, configInfo *par, struct grid *gp, molData *md\
  , imageInfo *img, double *lamtab, double *kaptab, const int nEntries){

  /*
This function constructs an image cube by following sets of rays (at least 1 per image pixel) through the model, solving the radiative transfer equations as appropriate for each ray. The ray locations within each pixel are chosen randomly within the pixel, but the number of rays per pixel is set equal to the number of projected model grid points falling within that pixel, down to a minimum equal to par->alias.

Note that the argument 'md', and the grid element '.mol', are only accessed for line images.
  */

  printf("\nRaytracing in progress...\n");

  const int minNumRaysForAverage=2;
  const int nsupsamppix = 3; // Number of pixels (in x,y) from origin in which to do cartesian supersampling
  const int supsamp = 10; // Number of rays per pixel (in x,y) for central pixel super-sampling

    double pixelSize, imgCentrePixels,minfreq,absDeltaFreq,xs[2],oneOnNumRays, ro;
    int totalNumImagePixels,ppi, ichan,lastChan,molI,lineI,i,j,di, xi, yi,id, index;
    double local_cmb,cmbFreq,scale;
    int cmbMolI,cmbLineI, ppx,ppy;

    pixelSize = img[im].distance*img[im].imgres;
    totalNumImagePixels = img[im].pxls*img[im].pxls;
    imgCentrePixels = img[im].pxls/2.0;


    if(img[im].doline){
    /* The user may have set img.trans/img.molI but not img.freq. If so, we calculate freq now.*/
      if(img[im].trans>-1)
        img[im].freq = md[img[im].molI].freq[img[im].trans];

      /* Fill in the missing one of the triplet nchan/velres/bandwidth.
      */
      if(img[im].bandwidth > 0 && img[im].velres > 0){
        img[im].nchan = (int)(img[im].bandwidth/(img[im].velres/CLIGHT*img[im].freq));

      }else if(img[im].bandwidth > 0 && img[im].nchan > 0){
        img[im].velres = img[im].bandwidth*CLIGHT/img[im].freq/img[im].nchan;

      }else{ /*(img[im].velres > 0 && img[im].nchan > 0 */
        img[im].bandwidth = img[im].nchan*img[im].velres/CLIGHT*img[im].freq;
      }
    } /* If not doline, we already have img.freq and nchan by now anyway. */

  /*
  We need to calculate or choose a single value of 'local' CMB flux, also single values (i.e. one of each per grid point) of dust and knu, all corresponding the the nominal image frequency. The sensible thing would seem to be to calculate them afresh for each new image; and for continuum images, this is what in fact has always been done. For line images however local_cmb and the dust/knu values were calculated for the frequency of each spectral line and stored respectively in the molData struct and the struct populations element of struct grid. These multiple values (of dust/knu at least) are required during the main solution kernel of LIME, so for line images at least they were kept until the present point, just so one from their number could be chosen. :-/
  At the present point in the code, for line images, instead of calculating the 'continuum' values of local_cmb/dust/knu, the algorithm chose the nearest 'line' frequency and calculates the required numbers from that. The intent is to preserve (for the present at least) the former numerical behaviour, while changing the way the information is parcelled out among the variables and structs. I.e. a dedicated 'continuum' pair of dust/knu values is now available for each grid point in addition to the array of spectral line values. This decoupling allows better management of memory and avoids the deceptive use of spectral-line variables for continuum use.
  */
   if(img[im].doline){
    if (img[im].trans>=0){
      cmbMolI  = img[im].molI;
      cmbLineI = img[im].trans;

    }else{ /* User didn't set trans. Find the nearest line to the image frequency. */
      minfreq = fabs(img[im].freq - md[0].freq[0]);;
      cmbMolI = 0;
      cmbLineI = 0;
      for(molI=0;molI<par->nSpecies;molI++){
        for(lineI=0;lineI<md[molI].nline;lineI++){
          if((molI==0 && lineI==0)) continue;

          absDeltaFreq = fabs(img[im].freq - md[molI].freq[lineI]);
          if(absDeltaFreq < minfreq){
            minfreq = absDeltaFreq;
            cmbMolI = molI;
            cmbLineI = lineI;
          }
        }
      }
    }
    cmbFreq = md[cmbMolI].freq[cmbLineI];

  }else{ /* continuum image */
    cmbFreq = img[im].freq;
  }

  local_cmb = planckfunc(cmbFreq,LOCAL_CMB_TEMP);
  calcGridContDustOpacity(par, cmbFreq, lamtab, kaptab, nEntries, gp); /* Reads gp attributes x, dens, and t and writes attributes cont.dust and cont.knu. */


  for(ppi=0;ppi<totalNumImagePixels;ppi++){
    for(ichan=0;ichan<img[im].nchan;ichan++){
      img[im].pixel[ppi].intense[ichan] = 0.0;
      img[im].pixel[ppi].tau[    ichan] = 0.0;
    }
  }

  for(ppi=0;ppi<totalNumImagePixels;ppi++)
    img[im].pixel[ppi].numRays = 0;

  double radiusarr[par->pIntensity], sorted_radius[par->pIntensity];

  for (id = 0;id < par->pIntensity;id++) {
    radiusarr[id] = sqrt(gp[id].x[0] * gp[id].x[0] + gp[id].x[1] * gp[id].x[1] + gp[id].x[2] * gp[id].x[2]);
    sorted_radius[id] = radiusarr[id]; //Note: not sorted yet
  }

  qsort(sorted_radius, par->pIntensity, sizeof(double), compare);

  struct radius_struct radius_struct;
  radius_struct.id = malloc(sizeof(int) * par->pIntensity);
  radius_struct.flux = malloc(sizeof(*radius_struct.flux) * par->pIntensity);
  radius_struct.radius = sorted_radius;
  double current;

  for (i = 0; i < par->pIntensity; i++) {
    current = radius_struct.radius[i];
    for (j = 0; j < par->pIntensity; j++) //TODO: More efficient algorithm than sequential search could be implemented
      if (current == radiusarr[j])
        radius_struct.id[i] = j; //holds sorted ids according to the time/radius of its corresponding gridpoint
  }

  for(i = 0; i < par->pIntensity; i++){
    radius_struct.flux[i].intense= malloc(sizeof(double) * img[im].nchan);
    radius_struct.flux[i].tau= malloc(sizeof(double) * img[im].nchan);

    for(ichan=0;ichan<img[im].nchan;ichan++){
      radius_struct.flux[i].intense[ichan] = 0.0;
      radius_struct.flux[i].tau[ichan] = 0.0;
    }
  }

  for(i = 0; i < par->pIntensity; i++){
    traceray(img,par,gp,md,radius_struct,im,i);
  }

  ppy = 0;
  for (ppi = 0; ppi<totalNumImagePixels; ppi++){
    ppx = ppi % img[im].pxls;
    if (ppx == 0 && ppi != 0)
      ppy += 1;

    xs[0] = (0.5 + ppx - imgCentrePixels) * pixelSize;
    xs[1] = (0.5 + ppy - imgCentrePixels) * pixelSize;
    ro = sqrt(xs[0]*xs[0] + xs[1]*xs[1]);

    for(i=0;i<par->pIntensity;i++) //TODO: use a more efficient searching algorithm, such as binary search (since its already sorted)
      if(radius_struct.radius[i]>ro){
        index =i;
        break;
      }
   // if(abs(xs[0])>pixelSize*nsupsamppix || abs(xs[1])>pixelSize*nsupsamppix){
      img[im].pixel[ppi].numRays++;

      if(index==0){
        for(ichan=0;ichan<img[im].nchan;ichan++){
          img[im].pixel[ppi].intense[ichan] = pow(10.0,radius_struct.flux[index].intense[ichan]);
          img[im].pixel[ppi].tau[ichan] = pow(10.0,radius_struct.flux[index].tau[ichan]);
        }
      }

      else if(index==par->pIntensity){
        index = par->pIntensity-1;
        for(ichan=0;ichan<img[im].nchan;ichan++){
          img[im].pixel[ppi].intense[ichan] = pow(10.0,radius_struct.flux[index].intense[ichan]);
          img[im].pixel[ppi].tau[ichan] = pow(10.0,radius_struct.flux[index].tau[ichan]);
        }
      }
      else{
        for(ichan=0;ichan<img[im].nchan;ichan++){
          img[im].pixel[ppi].intense[ichan] = pow(10.0,linear_interp(radius_struct.radius[index-1],radius_struct.radius[index],log10(radius_struct.flux[index-1].intense[ichan]),log10(radius_struct.flux[index].intense[ichan]),ro));
          img[im].pixel[ppi].tau[ichan] = pow(10.0,linear_interp(radius_struct.radius[index-1],radius_struct.radius[index],log10(radius_struct.flux[index-1].tau[ichan]),log10(radius_struct.flux[index].tau[ichan]),ro));
        }
      }
    //}
  }

  // scale = pixelSize/(double)supsamp;

  // for(j=1-(nsupsamppix*supsamp);j<(nsupsamppix*supsamp);j++){
  //   for(i=1-(nsupsamppix*supsamp);i<(nsupsamppix*supsamp);i++){
  //     xs[0] = j*scale;
  //     xs[1] = i*scale;
  //     ro = sqrt(xs[0]*xs[0] + xs[1]*xs[1]);

  //     xi = floor(xs[0]/pixelSize + imgCentrePixels);
  //     yi = floor(xs[1]/pixelSize + imgCentrePixels);
  //     ppi = yi * img[im].pxls + xi;
  //     img[im].pixel[ppi].numRays++;

  //     for(i=0;i<par->pIntensity;i++) //TODO: use a more efficient searching algorithm, such as binary search (since its already sorted)
  //       if(radius_struct.radius[i]>ro){
  //         index =i;
  //         break;
  //       }
  //     for(ichan=0;ichan<img[im].nchan;ichan++){
  //         img[im].pixel[ppi].intense[ichan] += pow(10.0,linear_interp(radius_struct.radius[index-1],radius_struct.radius[index],log10(radius_struct.flux[index-1].intense[ichan]),log10(radius_struct.flux[index].intense[ichan]),ro));
  //         img[im].pixel[ppi].tau[ichan] += pow(10.0,linear_interp(radius_struct.radius[index-1],radius_struct.radius[index],log10(radius_struct.flux[index-1].tau[ichan]),log10(radius_struct.flux[index].tau[ichan]),ro));
  //       }
  //   }
  // }

  for(ppi=0;ppi<totalNumImagePixels;ppi++){
    if(img[im].pixel[ppi].numRays >= minNumRaysForAverage){
      oneOnNumRays = 1.0/(double)img[im].pixel[ppi].numRays;
      for(ichan=0;ichan<img[im].nchan;ichan++){
        img[im].pixel[ppi].intense[ichan] *= oneOnNumRays;
        img[im].pixel[ppi].tau[    ichan] *= oneOnNumRays;
      }
    }
  }

  for(i = 0; i < par->pIntensity; i++){
    free(radius_struct.flux[i].intense);
    free(radius_struct.flux[i].tau);    
  }
  free(radius_struct.id);
  free(radius_struct.flux);


  if(par->polarization){ /* just add cmb to Stokes I, which is the first 'channel' */
    lastChan = 0;
  }else{
    lastChan = img[im].nchan;
  }

#ifdef FASTEXP
  for(ppi=0;ppi<totalNumImagePixels;ppi++){
    for(ichan=0;ichan<lastChan;ichan++){
      img[im].pixel[ppi].intense[ichan] += (FastExp(img[im].pixel[ppi].tau[ichan])-1.0)*local_cmb;
    }
  }
#else
  for(ppi=0;ppi<totalNumImagePixels;ppi++){
    for(ichan=0;ichan<lastChan;ichan++){
      img[im].pixel[ppi].intense[ichan] += (exp(   -img[im].pixel[ppi].tau[ichan])-1.0)*local_cmb;
    }
  }
#endif

}