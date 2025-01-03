/*
 *  sublime.h
 *  
 *  This file is part of SUBLIMED1D (SUBlimating gases in LIME, version D1D)
 *
 *  Copyright (C) 2006-2014 Christian Brinch
 *  Copyright (C) 2015-2017 The LIME development team
 *  Copyright (C) 2023 Martin Cordiner and Kristen Darnell (NASA GSFC)
 *
 */


#ifndef LIME_H
#define LIME_H

#ifdef IS_PYTHON
#include <Python.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <signal.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <fitsio.h>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_threads() 0
#define omp_get_thread_num() 0
#define omp_set_dynamic(int) 0
#endif

#include "dims.h"

#define VERSION "D1D-1.2"
#define DEFAULT_NTHREADS 1
#ifndef NTHREADS /* Value passed from the LIME script */
#define NTHREADS DEFAULT_NTHREADS
#endif

#include "constants.h"

#define HPIP            8.918502221e-27      /* HPLANCK*CLIGHT/4.0/PI/SPI	*/
#define HCKB            1.43877735           /* 100.*HPLANCK*CLIGHT/KBOLTZ	*/

/* Other constants */
#define NITERATIONS             16
#define MAX_RAYS_PER_POINT      10000
#define RAYS_PER_POINT          200
#define IMG_MIN_ALLOWED         1.0e-30
#define TOL                     1e-6
#define MAXITER                 50
#define maxBlendDeltaV          1.e4                   /* m/s */
#define N_RAN_PER_SEGMENT       3
#define FAST_EXP_MAX_TAYLOR     3
#define FAST_EXP_NUM_BITS       8
#define N_SMOOTH_ITERS          5                      /* number of smoothing iterations if using  par->samplingAlgorithm=0 */
#define TREE_POWER              2.0
#define ERF_TABLE_LIMIT         6.0                    /* For x>6 erf(x)-1<double precision machine epsilon, so no need to store the values for larger x. */
#define ERF_TABLE_SIZE          6145
#define BIN_WIDTH               (ERF_TABLE_LIMIT/(ERF_TABLE_SIZE-1.))
#define IBIN_WIDTH              (1./BIN_WIDTH)
#define N_VEL_SEG_PER_HALF      1
#define NUM_VEL_COEFFS          (1+2*N_VEL_SEG_PER_HALF) /* This is the number of velocity samples per edge (not including the grid vertices at each end of the edge). Currently this is elsewhere hard-wired at 3, the macro just being used in the file I/O modules. Note that we want an odd number of velocity samples per edge if we want to have the ability to do 2nd-order interpolation of velocity within Delaunay tetrahedra. */
#define MAX_NEG_OPT_DEPTH	30.0			/* 30 was the original value in LIME. */
#define NUM_RAN_DENS		100
#define RTOL 1.0e-10   /* CVODE scalar relative tolerance */
#define ATOL 1.0e-10   /* CVODE vector absolute tolerance components */
#define MINTOL 1.0e-15 /*Minimum allowed value for RTOL and ATOL*/

#define DEFAULT_BETA 1.042e-5 /*Standard water photolysis rate*/
#define DEFAULT_XNE 0.2 /*Standard scaling for the electron abundance*/
#define DEFAULT_TNUC 100.0 /*Default temperature for the nucleus*/
#define DEFAULT_PINTENSITY 300 /*Default number of points in radial grid points*/

/* Bit locations for the grid data-stage mask, that records the information which is present in the grid struct: */
#define DS_bit_x             0	/* id, x, sink */
#define DS_bit_neighbours    1	/* neigh, dir, ds, numNeigh */
#define DS_bit_velocity      2	/* vel */
#define DS_bit_density       3	/* dens */
#define DS_bit_abundance     4	/* abun, nmol */
#define DS_bit_turb_doppler  5	/* dopb */
#define DS_bit_temperatures  6	/* t */
#define DS_bit_magfield      7	/* B */
#define DS_bit_ACOEFF        8	/* a0, a1, a2, a3, a4 */
#define DS_bit_populations   9	/* mol */

#define DS_mask_x             1<<DS_bit_x
#define DS_mask_neighbours   (1<<DS_bit_neighbours   | DS_mask_x)
#define DS_mask_velocity     (1<<DS_bit_velocity     | DS_mask_x)
#define DS_mask_density      (1<<DS_bit_density      | DS_mask_x)
#define DS_mask_abundance    (1<<DS_bit_abundance    | DS_mask_x)
#define DS_mask_turb_doppler (1<<DS_bit_turb_doppler | DS_mask_x)
#define DS_mask_temperatures (1<<DS_bit_temperatures | DS_mask_x)
#define DS_mask_magfield     (1<<DS_bit_magfield     | DS_mask_x)
#define DS_mask_ACOEFF       (1<<DS_bit_ACOEFF       | DS_mask_neighbours | DS_mask_velocity)

#define DS_mask_1            DS_mask_x
#define DS_mask_2            DS_mask_neighbours
#define DS_mask_3            (DS_mask_2|DS_mask_density|DS_mask_temperatures)
#define DS_mask_4            (DS_mask_2|DS_mask_density|DS_mask_temperatures|DS_mask_abundance|DS_mask_turb_doppler|DS_mask_ACOEFF)
#define DS_mask_populations  (1<<DS_bit_populations | DS_mask_4)
#define DS_mask_5            DS_mask_populations
#define DS_mask_all          (DS_mask_populations | DS_mask_magfield)
#define DS_mask_all_but_mag  DS_mask_all & ~(1<<DS_bit_magfield)


#include "ufunc_types.h"
#include "collparts.h"
#include "inpars.h"
#include "defaults.h" /* includes lime_config.h */

struct cpData {
  double *down,*temp;
  int collPartId,ntemp,ntrans,*lcl,*lcu,densityIndex;
  char *name;
};

/* Previously declared in gridio.h. I moved it here so I could remove the unused file*/
struct linkType {
  unsigned int id, gis[2];
  double *vels;
};

struct molInfoType{
  char *molName;
  int nLevels, nLines;
};

struct gridInfoType{
  unsigned int nInternalPoints, nSinkPoints, nLinks, nNNIndices;
  unsigned short nDims, nSpecies, nDensities, nLinkVels;
  struct molInfoType *mols;
};

struct keywordType{
  int datatype; /* Codes given above. */
  char *keyname,*comment,*charValue;
  int intValue;
  float floatValue;
  double doubleValue;
  _Bool boolValue;
};

/* Molecular data: shared attributes */
typedef struct {
  int nlev,nline,npart;
  int *lal,*lau;
  double *aeinst,*freq,*beinstu,*beinstl,*eterm,*gstat,*gir;
  double *cmb,amass;
  struct cpData *part;
  char molName[80];
} molData;

/* Struct to store CKC data*/
typedef struct CKCdata{
	double ***Tedata; // [q][h][r] 
	double ***nedata; // [q][h][r] 
	int nr_values;
	double *r_values;
	int nQ_values;
	double *Q_values;
	int nH_values;
	double *H_values;
} CKCdata;

struct point {
  double x[DIM];
  double xn[DIM];
};

struct rates {
  int t_binlow;
  double interp_coeff;
};

struct continuumLine{
  double dust, knu;
};

struct populations {
  double *pops,*specNumDens;
  double dopb,binv,nmol,abun;
  struct rates *partner;
  struct continuumLine *cont;
};

/* Grid properties */
struct grid {
  int id;
  double x[DIM], vel[DIM], B[3]; /* B field only makes physical sense in 3 dimensions. */
  double *v1,*v2,*v3;
  double radius;
  int numNeigh;
  struct point *dir;
  struct grid **neigh;
  double *w;
  int sink;
  int nphot;
  int conv;
  double *dens,t[2],dopb_turb;
  double *ds;
  struct populations *mol;
  struct continuumLine cont;
};

struct spec {
  double *intense;
  double *tau;
  double stokes[3];
  int numRays;
};

/* Image information */
typedef struct {
  int doline;
  int nchan,trans,molI;
  struct spec *pixel;
  double velres;
  double imgres;
  double fwhm;
  int pxls;
  char *units;
  int *imgunits;
  int numunits;
  int psfShape;
  double psfWidth, binWidth;
  int rebinSpec;
  int nBins;
  double *psfKernel;
  int psfKernelN;
  double freq,bandwidth;
  char *filename;
  double source_vel;
  double theta,phi,incl,posang,azimuth;
  double distance;
  double rotMat[3][3];
  _Bool doInterpolateVels;
} imageInfo;

/* NOTE that it is assumed that vertx[i] is opposite the face that abuts with neigh[i] for all i.
*/ 
struct cell {
  struct grid *vertx[DIM+1];
  struct cell *neigh[DIM+1]; /* ==NULL flags an external face. */
  unsigned long id;
  double centre[DIM];
};

/* Some global variables */
extern _Bool fixRandomSeeds;

/* More functions */
int	run(inputPars, image*, const int);

_Bool	allBitsSet(const int flags, const int mask);
_Bool	anyBitSet(const int flags, const int mask);
_Bool	bitIsSet(const int flags, const int bitI);
_Bool	onlyBitsSet(const int flags, const int mask);

void	binpopsout(configInfo*, struct grid*, molData*);
void	calcDustData(configInfo*, double*, double*, const double, double*, const int, const double ts[], double*, double*);
void	calcExpTableEntries(const int, const int);
void	calcGridMolDensities(configInfo*, struct grid**);
void	calcGridMolDoppler(configInfo*, molData*, struct grid*);
void	calcGridMolSpecNumDens(configInfo*, molData*, struct grid*);
void	calcSourceFn(double, const configInfo*, double*, double*);
_Bool	charPtrIsNullOrEmpty(const char *inStr);
void	checkFirstLineMolDat(FILE *fp, char *moldatfile);
void	checkFgets(char *fgetsResult, char *message);
void	checkFscanf(const int fscanfResult, const int expectedNum, char *message);
void	checkFread(const size_t freadResult, const size_t expectedNum, char *message);
void	checkFwrite(const size_t fwriteResult, const size_t expectedNum, char *message);
void	checkUserDensWeights(configInfo*);
void	copyInparStr(const char*, char**);
void	delaunay(const int, struct grid*, const unsigned long, const _Bool, const _Bool, struct cell**, unsigned long*);
void	distCalc(configInfo*, struct grid*);
double	dotProduct3D(const double*, const double*);
double	FastExp(const float);
void	fillErfTable(void);
void	freeArrayOfStrings(char **arrayOfStrings, const int numStrings);
void	freeConfigInfo(configInfo*);
void	freeGrid(const unsigned int, const unsigned short, struct grid*);
void	freeImgInfo(const int, imageInfo*);
void	freeInputPars(inputPars *par);
void	freeMolData(const int, molData*);
void	freePopulation(const unsigned short, struct populations*);
void	freeSomeGridFields(const unsigned int, const unsigned short, struct grid*);
double	gaussline(const double, const double);
double	geterf(const double, const double);
void	getEdgeVelocities(configInfo *, struct grid *);
void	input(inputPars*, image*, char*);
double	interpolateKappa(const double, double*, double*, const int, gsl_spline*, gsl_interp_accel*);
int	levelPops(molData*, configInfo*, struct grid*, int*, double*, double*, const int, double*, int*);
void	mallocAndSetDefaultGrid(struct grid**, const size_t, const size_t);
void	mallocAndSetDefaultMolData(const int, molData**);
void	molInit(configInfo*, molData*);
double	planckfunc(const double, const double);
void	popsout(configInfo*, struct grid*, molData*);
void	processFitsError(int);
void	raytrace(int, configInfo*, struct grid*, molData*, imageInfo*, double*, double*, const int, double*);
void	readDustFile(char*, double**, double**, int*);
void	readMolData(configInfo *par, molData *md, int **allUniqueCollPartIds, int *numUniqueCollPartsFound);
void	readOrBuildGrid(configInfo*, struct grid**);
unsigned long reorderGrid(const unsigned long, struct grid*);
void	reportInfsAtOrigin(const int, const double*, const char*);
void	setCollPartsDefaults(struct cpData*);
int	setupAndWriteGrid(configInfo *par, struct grid *gp, molData *md, char *outFileName);
void	setUpDensityAux(configInfo*, int*, const int);
void	sigintHandler(int sigI);
int   compare(const void *a,const void *b);
double linear_interp(double x0, double x1, double y0, double y1, double value);
void	sourceFunc_line(const molData*, const double, const struct populations*, const int, double*, double*);
void	sourceFunc_cont(const struct continuumLine, double*, double*);
void	sourceFunc_pol(double*, const struct continuumLine, double (*rotMat)[3], double*, double*);
void	writeFitsAllUnits(const int, configInfo*, imageInfo*);
void	writeGridIfRequired(configInfo*, struct grid*, molData*, const int);
void	write_VTK_unstructured_Points(configInfo*, struct grid*);
struct CKCdata readCKCdata(CKCdata *st);
struct CKCdata readCKCfile(CKCdata *st, char *Tefilename, char *nefilename, double Q_values, double H_values);
double get_Telec (CKCdata *st, double Q, double rH, double radius);
double get_nelec (CKCdata *st, double Q, double rH, double radius);
void remove_spaces (char* restrict str_trimmed, const char* restrict str_untrimmed);
int indexOf(char *elm, char *arr[], int arr_cnt);


/* Curses functions */

void	bail_out(char*);
void	casaStyleProgressBar(const int, int);
void	collpartmesg(char*, int);
void	collpartmesg2(char*);
void	collpartmesg3(int, int);
void	error(char*);
void	goodnight(struct timeval);

void	greetings(void);
void	greetings_parallel(int);
void	printDone(int);
void	printMessage(char *);
void	progressbar(double, int);
void	progressbar2(configInfo*, int, int, double, double, double);
void	reportOutput(char*);
void	screenInfo(void);
void	warning(char*);

#ifdef FASTEXP
extern double EXP_TABLE_2D[128][10];
extern double EXP_TABLE_3D[256][2][10];
extern double oneOver_i[FAST_EXP_MAX_TAYLOR+1];
#endif

extern double ERF_TABLE[ERF_TABLE_SIZE];

#endif /* LIME_H */

