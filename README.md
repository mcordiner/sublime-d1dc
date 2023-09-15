# SUBLIMED1DC

SUBlimated cometary gases in LIME (Dynamical, 1D, compiled version)

This is SUBLIMED1D, a 1D radiative transfer code for outflowing cometary gases, by Martin Cordiner, Emmanuel Garcia-Berrios and Kristen Darnell (2023). 

SUBLIMED1DC is based on SUBLIMED, and the LIME (LIne Modeling Engine) version 1.9.3 by Christian Brinch (2006-2014) and the LIME development team (2015-2018). Main changes from the original LIME code include (1) switching from the static (GSL) solver to a time-dependent (CVODE) solver, which allows the dynamical nature of the cometary outflow to be properly simulated. Equations of statistical equilibrium in solver.c have been restructured as time-dependent differential equations, assuming constant outflow velocity. (2) Electron collision rates are added to the differential equations; analytic electron densities and temperatures are generated at runtime. (3) This model is strictly 1D, and the raytracing is now performed along a radial vector from the origin (in the plane of the sky), which is then interpolated onto the 2D image grid, supersampling the central pixels where the brightness can change rapidly. (4) For photon trapping, the escape probability approximation can invoked using par->useEP=1. (5) This version is compiled once, and the model parameters are read in at runtime using the "input.par" file; a few other parameters can be set in src/model.c before compiling.

Example input and output files are given in the example/ folder. 

To install, you will require the  gsl, qhull, cfitsio and cvode version 5 (SUNDIALS) libraries. If you have them installed outside the usual system lib folders, you will need to insert the relevant folder names in the Makefile.defs files, along with the correct qhull name (qhull|qhulstatic) and C compiler.

The code is then compiled by running the following shell script:

> ./compile_sublimed1d

This generates a binary executable called sublimed1dc.

To generate a reliable coma image, care needs to be taken to correctly set the par->radius parameter in model.c to capture all the expected emission (this will vary depending on the molecule, coma and viewing geometry). The channel spacing (velres) needs to be small enough (around 0.1 km/s or less) to properly sample the spectral line profile, even if the velocity information is later discarded. 

Experience has shown that a good model can be produced with par->pIntensity ~ 300 - 500 radial grid points. Qhull is used to generate the radial (1D) grid, weighted by the density. Radiation trapping effects tend to be very small, so useEP = 0 can be set for most models (apart from H2O), which allows the code to run much faster (in a matter of seconds), particularly for CH3OH. Faster raytracing is achieved for lower values of par->pIntensity.

If CVODE has an error but the input model appears physically reasonable, it can usually be fixed by adjusting RTOL and ATOL. Their default values are both 1.0E-9, but some CH3OH models require RTOL = ATOL = 1.0E-11, and some CO models require RTOL = ATOL = 1.0E-12. Some models may require RTOL to be increased to 1.0E-6, while ATOL can remain at 1.0E-9.


CKC_electrons branch: In this branch we can use the input parameter par->useCKC to control how the electron density and temperatures are determined. par->useCKC=0 uses the formalism of Zakharov et al. (2007). par->useCKC=1 interpolates the electron collision rates from data tables included in the program that were generated with the coma kinetic code developed by Dr. Cordiner and Dr. Charnley. par->useCKC=2 interpolates the electron collision rates from data tables provided by the user using par->CKCTeFile and par->CKCneFile.