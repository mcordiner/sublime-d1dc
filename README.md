# SUBLIME1D

SUBlimated Dynamical cometary gases in LIME

This is the first working version of SUBLIMED, a 3D Monte Carlo radiative transfer code for cometary comae, by Martin Cordiner and Miguel de Val Borro (2018). 

SUBLIMED is based on LIME (LIne Modeling Engine) version 1.9.3 by Christian Brinch (2006-2014) and the LIME development team (2015-2018). Main changes from the original LIME code include (1) the addition of electron collision rates to the static matrix in solver.c; analytic electron densities and temperatures are generated at runtime in solver.c; (2) the raytracing routine in raytrace.c has been altered to provide sufficient (evenly-weighted) sampling in the central image pixels to account for the strong, compact, central brightness peak of the coma.

Example input and output files are given in the example/ folder.

After installation of the required gsl and qhull libraries, recommend running the executable as follows:

> lime -f -n -p <#threads> model.c

To generate a reliable coma image, care needs to be taken to correctly set the par->radius parameter in model.c to capture all the expected emission (this will vary depending on the molecule, coma and viewing geometry). The channel spacing (velres) needs to be small enough (around 0.1 km/s or less) to properly sample the spectral line profile, even if the velocity information is later discarded. 

Experience has shown that a reasonably good model can be produced with only par->pIntensity = 1000 grid points (and par->sinkPoints = 500 sink points), with par->nSolveIters = 7. Such a coarsely-sampled model runs in only about 10 seconds on a 6-core workstation. For a more accurate model, use par->pIntensity = 5000, par->sinkPoints = 1000 and par->nSolveIters = 8 (to ensure convergence).
