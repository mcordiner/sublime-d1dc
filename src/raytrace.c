#include "lime.h"
#include "raythrucells.h"
#include <math.h>

struct flux{
  double *intense;
  double *tau;
};

struct molecule{
  int numLines;
  int *lines;
};

struct rayData{
  int *id;
  double *radius;
  struct molecule *mols;
  struct flux *flux;
};

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
traceray(imageInfo *img,configInfo *par,struct grid *gp,molData *md,struct rayData rayData,const int im, int index, double* rho_grid, int dz_grid_size, double* dz_grid, double* dz_vals, int* dz_indices, double* posneg){

// Note: this can be sped up considerably by separating the main loop into 2 parts: (1) generating the source function for each channel for each z point, and (2) doing the linear interpolation in a separate loop

  int ichan,stokesId,di,i,posn1,posn2,molI,lineI,lineID,zp_i;
  double zp,x[DIM],dx[DIM],vel[DIM],dz,dtau,col,r;
  double contJnu1,contAlpha1,contJnu2,contAlpha2,alpha,jnu,jnu1,alpha1,jnu2,alpha2,r1,r2,lineRedShift,vThisChan,deltav,logdz,vfac=0.0;
  double remnantSnu,expDTau,brightnessIncrement;

  for(di=0;di<DIM;di++){
    dx[di]= -img[im].rotMat[di][2]; /* This points towards the observer. */
  }

  x[0] = rho_grid[index];
  x[1] = 0.0;

  
  // Loop through the z-axis of the coma (from back to front), integrating the source function
  for (zp_i=0;zp_i<(2*dz_grid_size);zp_i++){
  
    x[2] = posneg[zp_i] * dz_grid[dz_indices[zp_i]];
    dz = dz_vals[dz_indices[zp_i]];
    r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    
    
    // Find the CVODE grid points that bracket our current radius, and if we are inside the model boundary, add to the integral
    if (r < rayData.radius[par->pIntensity - 1] && r > par->minScale){
      if (r > rayData.radius[0]){
       for (i = 1; i < par->pIntensity; i++){
         if (rayData.radius[i] > r){
             posn1 = rayData.id[i];
             posn2 = rayData.id[i-1];
             r1 = rayData.radius[i];
             r2 = rayData.radius[i-1];
             break;
         }
       }
      
      // Find the velocity vector   
      velocity(x[0],x[1],x[2],vel);
     
      /* Calculate first the continuum stuff because it is the same for all channels:*/
      contJnu1 = 0.0;
      contAlpha1 = 0.0; 
      sourceFunc_cont(gp[posn1].cont, &contJnu1, &contAlpha1);

      contJnu2 = 0.0;
      contAlpha2 = 0.0; 
      sourceFunc_cont(gp[posn2].cont, &contJnu2, &contAlpha2);

      for(ichan=0;ichan<img[im].nchan;ichan++){
        jnu1 = contJnu1;
        alpha1 = contAlpha1;
        jnu2 = contJnu2;
        alpha2 = contAlpha2;

        vThisChan = (ichan-(img[im].nchan-1)*0.5)*img[im].velres; /* Consistent with the WCS definition in writefits(). */
        if(img[im].doline){
          for(molI=0;molI<par->nSpecies;molI++){
            for(lineID=0;lineID<rayData.mols[molI].numLines;lineID++){
                lineI = rayData.mols[molI].lines[lineID];
                /* Calculate the red shift of the transition wrt to the frequency specified for the image.
                */
                if(img[im].trans > -1){
                  lineRedShift=(md[img[im].molI].freq[img[im].trans]-md[molI].freq[lineI])/md[img[im].molI].freq[img[im].trans]*CLIGHT;
                } else {
                  lineRedShift=(img[im].freq-md[molI].freq[lineI])/img[im].freq*CLIGHT;
                }
                deltav = vThisChan - img[im].source_vel - lineRedShift;

                /* Calculating the source function for the nearest two radial points, which will later be used to intepolate the jnu and alpha at the current radial point
                */

                //Calculating source function and velocity term for 1st radial point
                vfac = gaussline(deltav-dotProduct3D(dx,vel),gp[posn1].mol[molI].binv);                  
                sourceFunc_line(&md[molI],vfac,&(gp[posn1].mol[molI]),lineI,&jnu1,&alpha1);

                //Calculating source function and velocity term for 2nd radial point
                vfac = gaussline(deltav-dotProduct3D(dx,vel),gp[posn2].mol[molI].binv);       
                sourceFunc_line(&md[molI],vfac,&(gp[posn2].mol[molI]),lineI,&jnu2,&alpha2);

            }//end for lineI
          }//end for molI
        }//end if doline

        alpha = linear_interp(r1, r2, alpha1, alpha2, r);
        jnu = linear_interp(r1, r2, jnu1, jnu2, r);

        dtau=alpha*dz;
        calcSourceFn(dtau, par, &remnantSnu, &expDTau);
        remnantSnu *= jnu*dz;

#ifdef FASTEXP
        brightnessIncrement = FastExp(rayData.flux[index].tau[ichan])*remnantSnu;
#else
        brightnessIncrement =    exp(-rayData.flux[index].tau[ichan])*remnantSnu;
#endif

        rayData.flux[index].intense[ichan] += brightnessIncrement;
        rayData.flux[index].tau[ichan] += dtau;
      
      }//end for ichan
     }//check for r<radius[0]
   }//Move to next z point

  }//Loop over z points

}

/*....................................................................*/
void rhoGrid2image(int ppi,configInfo *par, double *rho_grid, struct rayData rayData, imageInfo *img, double pixelSize, double imgCentrePixels, const int im){
/*Take the raytraced rho vector and interpolate it into the 2D image grid at pixel number ppi*/
    double xs[2], ro;
    int ppx, ppy, ichan, i, index;
    
    ppx = ppi % img[im].pxls;
    ppy = (int)(ppi/img[im].pxls);

    xs[0] = (0.5 + ppx - imgCentrePixels) * pixelSize;
    xs[1] = (0.5 + ppy - imgCentrePixels) * pixelSize;
    ro = sqrt(xs[0]*xs[0] + xs[1]*xs[1]);
    
    if(ro>=rho_grid[par->pIntensity-1]){
    // Pixel is outside model domain so set to zero
      for(ichan=0;ichan<img[im].nchan;ichan++){
          img[im].pixel[ppi].intense[ichan] = 0.0;
          img[im].pixel[ppi].tau[ichan] = 0.0;          
        }
    }else{
    
    for(i=0;i<par->pIntensity;i++) //TODO: use a more efficient searching algorithm, such as binary search (since its already sorted)
      if(rho_grid[i]>ro){
        index =i;
        break;
      }

      if(index==0){
        for(ichan=0;ichan<img[im].nchan;ichan++){
          img[im].pixel[ppi].intense[ichan] = rayData.flux[index].intense[ichan];
          img[im].pixel[ppi].tau[ichan] = rayData.flux[index].tau[ichan];          
        }
      }

      else if(index==par->pIntensity){
        index = par->pIntensity-1;
        for(ichan=0;ichan<img[im].nchan;ichan++){
          img[im].pixel[ppi].intense[ichan] = rayData.flux[index].intense[ichan];
          img[im].pixel[ppi].tau[ichan] = rayData.flux[index].tau[ichan];
        }
      }
      else{
        for(ichan=0;ichan<img[im].nchan;ichan++){
          img[im].pixel[ppi].intense[ichan] = pow(10.0,linear_interp(rho_grid[index-1],rho_grid[index],log10(rayData.flux[index-1].intense[ichan]),log10(rayData.flux[index].intense[ichan]),ro));
          img[im].pixel[ppi].tau[ichan] = linear_interp(rho_grid[index-1],rho_grid[index],rayData.flux[index-1].tau[ichan],rayData.flux[index].tau[ichan],ro);
        }
      }
    }
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

    double pixelSize, imgCentrePixels,minfreq,absDeltaFreq,xs[2],rho_grid[par->pIntensity],radiusarr[par->pIntensity], sorted_radius[par->pIntensity], ro;
    int totalNumImagePixels,ppi, ichan,lastChan,molI,lineI,i,j,k,di, xi, yi,id, index, pixoff,pixoff2,pixshiftx,pixshifty, nsupsamppix,numRays;
    double local_cmb,cmbFreq,scale,shift,offset,logdz,*zp_grid = NULL,*dz_grid = NULL,*dz_vals = NULL,*posneg = NULL;
    int cmbMolI,cmbLineI, ppx,ppy,*dz_indices = NULL;

    // Set up the z integration grid. Length of the z grid is double this:
    int dz_grid_size = 50;
  
    // Image parameters
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

  for (id = 0;id < par->pIntensity;id++) {
    radiusarr[id] = sqrt(gp[id].x[0] * gp[id].x[0] + gp[id].x[1] * gp[id].x[1] + gp[id].x[2] * gp[id].x[2]);
    sorted_radius[id] = radiusarr[id]; //Note: not sorted yet
  }

  qsort(sorted_radius, par->pIntensity, sizeof(double), compare);

  struct rayData rayData;
  rayData.id = malloc(sizeof(int) * par->pIntensity);
  rayData.flux = malloc(sizeof(*rayData.flux) * par->pIntensity);
  rayData.radius = sorted_radius;
  double current;
  
  offset =  rayData.radius[0];
  
  for (i = 0; i < par->pIntensity; i++) {
    rho_grid[i] = rayData.radius[i] - offset;
  }

  for (i = 0; i < par->pIntensity; i++) {
    current = rayData.radius[i];
    for (j = 0; j < par->pIntensity; j++) //TODO: More efficient algorithm than sequential search could be implemented
      if (current == radiusarr[j])
        rayData.id[i] = j; //holds sorted ids according to the time/radius of its corresponding gridpoint
  }

  for(i = 0; i < par->pIntensity; i++){
    rayData.flux[i].intense= malloc(sizeof(double) * img[im].nchan);
    rayData.flux[i].tau= malloc(sizeof(double) * img[im].nchan);

    for(ichan=0;ichan<img[im].nchan;ichan++){
      rayData.flux[i].intense[ichan] = 0.0;
      rayData.flux[i].tau[ichan] = 0.0;
    }
  }

  rayData.mols = malloc(sizeof(*rayData.mols) * par->nSpecies);
  int numLines= 0;

  //For each molecule,count how many lines fall within the spectral range of the image. This way, we know how much space to allocate for each rayData.mols[molI].lines
  //We do these steps now so that we don't need to determine which spectral lines contribute to the image repeteadly inside traceray(), which is time consuming
  for(molI=0;molI<par->nSpecies;molI++){
    for(lineI=0;lineI<md[molI].nline;lineI++){
      if(md[molI].freq[lineI] > img[im].freq-img[im].bandwidth*0.5\
        && md[molI].freq[lineI] < img[im].freq+img[im].bandwidth*0.5){
        numLines++;
      }
    }
    rayData.mols[molI].lines = malloc(sizeof(*rayData.mols[molI].lines)*numLines);
    rayData.mols[molI].numLines = numLines;
    numLines = 0;
  }

  //Now we store the index of the lines that contribute to the image.
  for(molI=0;molI<par->nSpecies;molI++){
    index = 0;
    for(lineI=0;lineI<md[molI].nline;lineI++){
      if(md[molI].freq[lineI] > img[im].freq-img[im].bandwidth*0.5\
        && md[molI].freq[lineI] < img[im].freq+img[im].bandwidth*0.5){
        rayData.mols[molI].lines[index] = lineI;
        index++;
      }
    }
  }
  
  // Set up the z integration grid
  dz_grid = malloc(sizeof(*dz_grid) * dz_grid_size);
  dz_vals = malloc(sizeof(*dz_vals) * dz_grid_size);
  dz_indices = malloc(sizeof(*dz_indices) * dz_grid_size * 2);
  posneg = malloc(sizeof(*posneg) * dz_grid_size * 2);
  zp_grid = malloc(sizeof(*zp_grid) * dz_grid_size + 1);
  
   
  //Creating zp_grid
  logdz=log10(par->radius)/dz_grid_size;
  for (i=0;i<=dz_grid_size;i++){
  zp_grid[i]=pow(10,i*logdz)-1.0;
  }
  zp_grid[dz_grid_size]=par->radius;
  
  //Halfway points and associated dz values
  for (i=0;i<dz_grid_size;i++){
  dz_grid[i]=(zp_grid[i]+zp_grid[i+1])/2.0;
  dz_vals[i]=zp_grid[i+1]-zp_grid[i];
  }

  //Generate the indices and posneg terms (for front vs. back of coma)
  j=0;
  for (i=(dz_grid_size-1);i>=0;i--){
     dz_indices[j]=i;
     posneg[j]=1.0;
     j++;
  }
    for (i=0;i<dz_grid_size;i++){
     dz_indices[j]=i;
     posneg[j]=-1.0;
     j++;
  }

  //printf("Calling traceray...\n");
  
  // Do the raytracing as a function of rho (linear vector from the origin), to later be interpolated onto the image grid
  // Parallel loop over rho grid points
  omp_set_num_threads(par->nThreads);
  #pragma omp parallel for schedule (dynamic)
  for(i = 0; i < par->pIntensity; i++){
    traceray(img,par,gp,md,rayData,im,i,rho_grid,dz_grid_size,dz_grid,dz_vals,dz_indices,posneg);
  }

  free(dz_grid);
  free(dz_vals);
  free(zp_grid);
  free(dz_indices);
  
  //printf("Interpolating to image grid...\n");

  // Parallel loop over image pixels
  #pragma omp parallel for schedule (dynamic)
  for (ppi = 0; ppi<totalNumImagePixels; ppi++){
     rhoGrid2image(ppi,par,rho_grid,rayData,img,pixelSize,imgCentrePixels,im);
  }

  //printf("Supersampling the central pixels...\n");
  
  const int supsamp = 20; // Number of rays per pixel (supsamp * supsamp in x,y plane)
  scale = pixelSize/((double)supsamp);
  
  if(img[im].pxls % 2 != 0){
  // If there is an odd number of image pixels, supersample the innermost nsupsamppix x nsupsamppix region:
      nsupsamppix = 5; 
      shift = (pixelSize/2.0) + (scale/2.0);
      pixoff = 1;
      pixoff2 = 0;
  }else{
  // If there is an even number of image pixels, supersample the innermost nsupsamppix x nsupsamppix region:
      nsupsamppix = 4;
      shift = (scale/2.0);
      pixoff = 0;
      pixoff2 = 1;
  }
  
  for(pixshiftx=(pixoff-nsupsamppix)/2;pixshiftx<=(nsupsamppix-pixoff-pixoff2)/2;pixshiftx++){
  // Parallel loop over y axis of supersampled region. At some point we should extend this treatment to all pixels
    #pragma omp parallel for schedule (dynamic) default(shared) private(numRays,xi,yi,ppi,ichan,i,j,xs,ro,index,k)
    for(pixshifty=(pixoff-nsupsamppix)/2;pixshifty<=(nsupsamppix-pixoff-pixoff2)/2;pixshifty++){
      
      numRays=0;
      
      // Initially set supersampled pixels to zero
      xi=pixshiftx+(img[im].pxls-pixoff)/2;
      yi=pixshifty+(img[im].pxls-pixoff)/2;
      ppi = yi * img[im].pxls + xi;
      for(ichan=0;ichan<img[im].nchan;ichan++){
          img[im].pixel[ppi].intense[ichan] = 0.0;
          img[im].pixel[ppi].tau[ichan] = 0.0;
      }
      
      for(j=1;j<=supsamp;j++){
        for(i=1;i<=supsamp;i++){
          xs[0] = (j*scale) - shift + (pixshiftx*pixelSize);
          xs[1] = (i*scale) - shift + (pixshifty*pixelSize);
          ro = sqrt(xs[0]*xs[0] + xs[1]*xs[1]);
          
          if(ro>=par->radius){
            // Outside the model boundary so do not add any flux
            numRays++;
            continue;
          
          }else{
          
          xi = round(xs[0]/pixelSize + imgCentrePixels - 0.5);
          yi = round(xs[1]/pixelSize + imgCentrePixels - 0.5);
          ppi = yi * img[im].pxls + xi;
          numRays++;

          for(k=0;k<par->pIntensity;k++)
            if(rho_grid[k]>ro){
              index =k;
              break;
            }
          for(ichan=0;ichan<img[im].nchan;ichan++){
              img[im].pixel[ppi].intense[ichan] += pow(10.0,linear_interp(rho_grid[index-1],rho_grid[index],log10(rayData.flux[index-1].intense[ichan]),log10(rayData.flux[index].intense[ichan]),ro));
              img[im].pixel[ppi].tau[ichan] += linear_interp(rho_grid[index-1],rho_grid[index],rayData.flux[index-1].tau[ichan],rayData.flux[index].tau[ichan],ro);
            }
          }
        }
      }
    
     for(ichan=0;ichan<img[im].nchan;ichan++){
       img[im].pixel[ppi].intense[ichan] /= numRays;
       img[im].pixel[ppi].tau[ichan] /= numRays;
     }
    }
  }
  

  for(i = 0; i < par->pIntensity; i++){
    free(rayData.flux[i].intense);
    free(rayData.flux[i].tau);    
  }
  for(molI=0;molI<par->nSpecies;molI++){
    free(rayData.mols[molI].lines);
  }
  free(rayData.id);
  free(rayData.flux);
  free(rayData.mols);

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