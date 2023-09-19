# SUBLIME-D1DC

SUBlimated cometary gases in LIME (Dynamical, 1D, Compiled version)

This is SUBLIME-D1DC, a 1D radiative transfer code for outflowing cometary gases, by Martin Cordiner, Emmanuel Garcia-Berrios and Kristen Darnell (2023). 

SUBLIME-D1DC is based on the [SUBLIME code](https://ui.adsabs.harvard.edu/abs/2022ApJ...929...38C/abstract) by Martin Cordiner (2022), and the [LIME (LIne Modeling Engine)](https://github.com/lime-rt/lime) by Christian Brinch (2006-2014) and the LIME development team (2015-2018). Main changes from the original LIME code include (1) switching from the static (GSL) solver to a time-dependent (CVODE) solver, which allows the dynamical nature of the cometary outflow to be properly simulated. Equations of statistical equilibrium in solver.c have been restructured as time-dependent differential equations, assuming constant outflow velocity. (2) Electron collision rates are added to the differential equations; analytic electron densities and temperatures are generated at runtime. (3) This model is strictly 1D, and the raytracing is now performed along a radial vector from the origin (in the plane of the sky), which is then interpolated onto the 2D image grid, supersampling the central pixels where the brightness can change rapidly. (4) For photon trapping, the escape probability approximation can invoked using par->useEP=1. (5) This version is compiled once, and the model parameters are read in at runtime using the "input.par" file; a few other parameters can be set in src/model.c before compiling (see also sublime.h for the values of physical constants and other default parameters).

Example input and output files are given in the example/ folder. 

To install, you will require the  gsl, cfitsio and cvode version 5 (SUNDIALS) libraries. If you have them installed outside the usual system lib folders, you will need to insert the relevant folder names in the Makefile.defs files, along with the correct path to your C compiler.

The code is then compiled by running the following shell script:

> ./compile_sublimed1d

This generates a binary executable called sublimed1dc, which reads parameters from the input.par file (located in the current working directory). 

To generate a reliable coma image, care needs to be taken to correctly set the par->radius parameter in model.c to capture all the expected emission (this will vary depending on the molecule, coma and viewing geometry). The channel spacing (velres) needs to be small enough (around 100 m/s or less) to properly sample the spectral line profile, even if the velocity information is later discarded. 

The radial grid is generated pseudo-randomly, with a density of points proportional to the gas density. Experience has shown that a good model can be produced with npts ~ 300 - 500 (300 is the default).  Radiation trapping effects tend to be very small for molecules other than water, so useEP = 0 should be used for most molecules, which allows the code to run faster. Faster raytracing is achieved for lower values of npts, nchan and pxls.

If CVODE has an error but the input model appears physically reasonable, it can usually be fixed by adjusting RTOL and ATOL in sublime.h (the relative and absolute tolerances for parameter errors). Their default values are both 1.0E-10, but some CH3OH models require RTOL = ATOL = 1.0E-11, and some CO models require RTOL = ATOL = 1.0E-12. In the event that CVODE fails, the code will automatically try reducing RTOL and ATOL until success is achieved. Some models may require RTOL to be increased to 1.0E-6, while ATOL can remain at 1.0E-9.

The optional input parameter par->useCKCdata controls how the electron densities and temperatures are determined. par->useCKCdata=0 uses the formalism of Zakharov et al. (2007, A&A, 473, 303). par->useCKCdata=1 interpolates the electron collision rates from data tables included in the program that were generated with the coma kinetic code developed by Cordiner & Charnley (2021, MNRAS, 504, 5401). par->useCKC=2 interpolates the electron collision rates from data tables provided by the user using par->CKCTeFile and par->CKCneFile.

If you use this code in a published work, please reference Cordiner, M. A., Coulson, I. M., Garcia-Berrios, E. et al. 2022, Astrophysical Journal, Volume 929, id.38.
