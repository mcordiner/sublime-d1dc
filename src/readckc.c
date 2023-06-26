/*
 *  readckc.c
 *
 *  This file is part of SUBLIMED1D (SUBlimating gases in LIME, version D1D)
 *
 *  Copyright (C) 2023 Martin Cordiner and Kristen Darnell (NASA GSFC)
 *
 */




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lime.h"

/* FUNCTION TO READ CKC DATA: LOOKUP TABLES */
struct CKCdata readCKCdata(CKCdata *st){

	char *ParentDircPath = getenv("SUBLIMEDATA");
	//char *ParentDircPath = "/Users/kdarnell/code/sublimed1d/data/";
	char FILEPATH[250];
	
	int nfiles=35;
	char *Tefilename[35] = { 
		"ckc_data/Te_Q1e+27_H0.5.dat", 
		"ckc_data/Te_Q1e+27_H1.0.dat", 
		"ckc_data/Te_Q1e+27_H1.5.dat", 
		"ckc_data/Te_Q1e+27_H2.0.dat", 
		"ckc_data/Te_Q1e+27_H2.5.dat", 
		"ckc_data/Te_Q5e+27_H0.5.dat", 
		"ckc_data/Te_Q5e+27_H1.0.dat", 
		"ckc_data/Te_Q5e+27_H1.5.dat", 
		"ckc_data/Te_Q5e+27_H2.0.dat", 
		"ckc_data/Te_Q5e+27_H2.5.dat", 
		"ckc_data/Te_Q1e+28_H0.5.dat", 
		"ckc_data/Te_Q1e+28_H1.0.dat", 
		"ckc_data/Te_Q1e+28_H1.5.dat", 
		"ckc_data/Te_Q1e+28_H2.0.dat", 
		"ckc_data/Te_Q1e+28_H2.5.dat", 
		"ckc_data/Te_Q5e+28_H0.5.dat", 
		"ckc_data/Te_Q5e+28_H1.0.dat", 
		"ckc_data/Te_Q5e+28_H1.5.dat", 
		"ckc_data/Te_Q5e+28_H2.0.dat", 
		"ckc_data/Te_Q5e+28_H2.5.dat", 
		"ckc_data/Te_Q1e+29_H0.5.dat", 
		"ckc_data/Te_Q1e+29_H1.0.dat", 
		"ckc_data/Te_Q1e+29_H1.5.dat", 
		"ckc_data/Te_Q1e+29_H2.0.dat", 
		"ckc_data/Te_Q1e+29_H2.5.dat", 
		"ckc_data/Te_Q5e+29_H0.5.dat", 
		"ckc_data/Te_Q5e+29_H1.0.dat", 
		"ckc_data/Te_Q5e+29_H1.5.dat", 
		"ckc_data/Te_Q5e+29_H2.0.dat", 
		"ckc_data/Te_Q5e+29_H2.5.dat", 
		"ckc_data/Te_Q1e+30_H0.5.dat", 
		"ckc_data/Te_Q1e+30_H1.0.dat", 
		"ckc_data/Te_Q1e+30_H1.5.dat", 
		"ckc_data/Te_Q1e+30_H2.0.dat", 
		"ckc_data/Te_Q1e+30_H2.5.dat" 
		};
	char *nefilename[35] = { 
		"ckc_data/ne_Q1e+27_H0.5.dat", 
		"ckc_data/ne_Q1e+27_H1.0.dat", 
		"ckc_data/ne_Q1e+27_H1.5.dat", 
		"ckc_data/ne_Q1e+27_H2.0.dat", 
		"ckc_data/ne_Q1e+27_H2.5.dat", 
		"ckc_data/ne_Q5e+27_H0.5.dat", 
		"ckc_data/ne_Q5e+27_H1.0.dat", 
		"ckc_data/ne_Q5e+27_H1.5.dat", 
		"ckc_data/ne_Q5e+27_H2.0.dat", 
		"ckc_data/ne_Q5e+27_H2.5.dat", 
		"ckc_data/ne_Q1e+28_H0.5.dat", 
		"ckc_data/ne_Q1e+28_H1.0.dat", 
		"ckc_data/ne_Q1e+28_H1.5.dat", 
		"ckc_data/ne_Q1e+28_H2.0.dat", 
		"ckc_data/ne_Q1e+28_H2.5.dat", 
		"ckc_data/ne_Q5e+28_H0.5.dat", 
		"ckc_data/ne_Q5e+28_H1.0.dat", 
		"ckc_data/ne_Q5e+28_H1.5.dat", 
		"ckc_data/ne_Q5e+28_H2.0.dat", 
		"ckc_data/ne_Q5e+28_H2.5.dat", 
		"ckc_data/ne_Q1e+29_H0.5.dat", 
		"ckc_data/ne_Q1e+29_H1.0.dat", 
		"ckc_data/ne_Q1e+29_H1.5.dat", 
		"ckc_data/ne_Q1e+29_H2.0.dat", 
		"ckc_data/ne_Q1e+29_H2.5.dat", 
		"ckc_data/ne_Q5e+29_H0.5.dat", 
		"ckc_data/ne_Q5e+29_H1.0.dat", 
		"ckc_data/ne_Q5e+29_H1.5.dat", 
		"ckc_data/ne_Q5e+29_H2.0.dat", 
		"ckc_data/ne_Q5e+29_H2.5.dat", 
		"ckc_data/ne_Q1e+30_H0.5.dat", 
		"ckc_data/ne_Q1e+30_H1.0.dat", 
		"ckc_data/ne_Q1e+30_H1.5.dat", 
		"ckc_data/ne_Q1e+30_H2.0.dat", 
		"ckc_data/ne_Q1e+30_H2.5.dat" 
		};
	int nQ_values=7;
	double Q_values[7] = { 
		1.00e+27, 
		5.00e+27, 
		1.00e+28, 
		5.00e+28, 
		1.00e+29, 
		5.00e+29, 
		1.00e+30 
		};
	int nH_values=5;
	double H_values[5] = { 
		5.00e-01, 
		1.00e+00, 
		1.50e+00, 
		2.00e+00, 
		2.50e+00 
		};
	
	int nr_values;
	double *r_values;
	
   	double ***Tedata; // [q][h][r] 
	double ***nedata; // [q][h][r] 

   	int nrows=-1, ncolumns=1;
   	int headerrows=3; // 4-1 
   
   	int i, j, k, ir=0, iq=0, ih=0; //counters
   	char c;  // To store a character read from file
   
   	FILE * fp[nfiles];
   		
   	ncolumns = nfiles + 1;
  	
	// Counting the number of rowss
	strcpy(FILEPATH, ParentDircPath);
	strcat(FILEPATH, "/");
	strcat(FILEPATH, Tefilename[0]);
	
	fp[0] = fopen (FILEPATH, "r+");
	
	for (c = getc(fp[0]); c != EOF; c = getc(fp[0]))
    	if (c == '\n') // Increment count if this character is newline
       		nrows = nrows + 1;
    
    nr_values = nrows - headerrows -2; // there were two extra rows being added to the array so I decreased nrows by 2
    
	rewind(fp[0]); //Returns to start of file
	fclose(fp[0]);
	
	// allocate the memory for the array, code from: http://www.eskimo.com/~scs/cclass/int/sx9b.html
	r_values = malloc(nr_values * sizeof(double));
	
	Tedata = malloc(nQ_values * sizeof(double));
	nedata = malloc(nQ_values * sizeof(double));
	for(iq = 0; iq < nQ_values; iq++) {
		Tedata[iq] = malloc(nH_values * sizeof(double));
		nedata[iq] = malloc(nH_values * sizeof(double));
		for(ih = 0; ih < nH_values; ih++) {
			Tedata[iq][ih] = malloc(nr_values * sizeof(double));
			nedata[iq][ih] = malloc(nr_values * sizeof(double));
		}
	}
	
	
	// Reading the Te files
	ir = 0;
	iq = 0;
	ih = 0;
	
	for(int i=0; i < nfiles; i++) {
		strcpy(FILEPATH, ParentDircPath);
		strcat(FILEPATH, "/");
		strcat(FILEPATH, Tefilename[i]);
	
		fp[i] = fopen (FILEPATH, "r+");

		//Skips the header rows
		for(j = 0; j < headerrows; j++){
		    fscanf(fp[i], "%*[^\n]\n"); 
		}
	
		// Reads the data into the 2D array line by line
		for(ir = 0; ir < nr_values; ir++){
    		fscanf(fp[i], "%le %le", &r_values[ir], &Tedata[iq][ih][ir]);
  		}
 		
  		// Incriments and resets q and h indicies
  		ih++;
  		if (ih >= nH_values) {
  			ih = 0;
  			iq++;
  		} 		
  		if (iq >= nQ_values) iq=0;
  		
	  	fclose(fp[i]);
	  } 
	
	// Reading ne files
	ir = 0;
	iq = 0;
	ih = 0;
	
	for(int i=0; i < nfiles; i++) {
		strcpy(FILEPATH, ParentDircPath);
		strcat(FILEPATH, "/");
		strcat(FILEPATH, nefilename[i]);
		
		fp[i] = fopen (FILEPATH, "r+");

		//Skips the header rows
		for(j = 0; j < headerrows; j++){
		    fscanf(fp[i], "%*[^\n]\n"); 
		}
	
		// Reads the data into the 2D array line by line
		for(ir = 0; ir < nr_values; ir++){
    		fscanf(fp[i], "%le %le", &r_values[ir], &nedata[iq][ih][ir]);
  		}
  		
  		// Incriments and resets q and h indicies
  		ih++;
  		if (ih >= nH_values) {
  			ih = 0;
  			iq++;
  		} 		
  		if (iq >= nQ_values) iq=0;

	  	fclose(fp[i]);
	  }   
    
    //Saving values into the structure to be retuned
    
    st->Tedata = Tedata;
    st->nedata = nedata;
	st->nr_values = nr_values;
	st->r_values = r_values; 
	st->nQ_values = nQ_values;
	st->Q_values = Q_values;
	st->nH_values = nH_values;
	st->H_values = H_values;
	    
	return *st;
}

/* FUNCTION TO READ CKC DATA: SINGLE TABLE */
/* TODO: Remove unnecessary loops */
struct CKCdata readCKCfile(CKCdata *st, char *Tefilename, char *nefilename, double Q_values, double H_values){

	//char *ParentDircPath = getenv("SUBLIMEDATA");
	//char *ParentDircPath = "/Users/kdarnell/code/sublimed1d/data/";
	char FILEPATH[250];
	
	int nfiles=1;
	int nQ_values=1;
	int nH_values=1;
	
	int nr_values;
	double *r_values;
	
   	double ***Tedata; // [q][h][r] 
	double ***nedata; // [q][h][r] 

   	int nrows=-1, ncolumns=1;
   	int headerrows=3; // 4-1 
   
   	int i, j, k, ir=0, iq=0, ih=0; //counters
   	char c;  // To store a character read from file
   
   	FILE * fp[nfiles];
   		
   	ncolumns = nfiles + 1;
   	  	
	// Counting the number of rowss
	//strcpy(FILEPATH, ParentDircPath);
	//strcat(FILEPATH, "/");
	//strcat(FILEPATH, Tefilename);
	strcpy(FILEPATH, Tefilename);
	
	fp[0] = fopen (FILEPATH, "r+");
	
	for (c = getc(fp[0]); c != EOF; c = getc(fp[0]))
    	if (c == '\n') // Increment count if this character is newline
       		nrows = nrows + 1;
    
    nr_values = nrows - headerrows -2; // there were two extra rows being added to the array so I decreased nrows by 2
    //printf("nr_values = %d \n", nr_values);
    
	rewind(fp[0]); //Returns to start of file
	fclose(fp[0]);
	
	// allocate the memory for the array, code from: http://www.eskimo.com/~scs/cclass/int/sx9b.html
	r_values = malloc(nr_values * sizeof(double));
	
	Tedata = malloc(nQ_values * sizeof(double));
	nedata = malloc(nQ_values * sizeof(double));
	for(iq = 0; iq < nQ_values; iq++) {
		Tedata[iq] = malloc(nH_values * sizeof(double));
		nedata[iq] = malloc(nH_values * sizeof(double));
		for(ih = 0; ih < nH_values; ih++) {
			Tedata[iq][ih] = malloc(nr_values * sizeof(double));
			nedata[iq][ih] = malloc(nr_values * sizeof(double));
		}
	}
	
	
	// Reading the Te files
	ir = 0;
	iq = 0;
	ih = 0;
	
	for(int i=0; i < nfiles; i++) {
		//strcpy(FILEPATH, ParentDircPath);
		//strcat(FILEPATH, "/");
		//strcat(FILEPATH, Tefilename);
		strcpy(FILEPATH, Tefilename);
	
		fp[i] = fopen (FILEPATH, "r+");

		//Skips the header rows
		for(j = 0; j < headerrows; j++){
		    fscanf(fp[i], "%*[^\n]\n"); 
		}
	
		// Reads the data into the 2D array line by line
		for(ir = 0; ir < nr_values; ir++){
    		fscanf(fp[i], "%le %le", &r_values[ir], &Tedata[iq][ih][ir]);
  		}
 		
  		// Incriments and resets q and h indicies
  		ih++;
  		if (ih >= nH_values) {
  			ih = 0;
  			iq++;
  		} 		
  		if (iq >= nQ_values) iq=0;

	  	fclose(fp[i]);
	  } 
	
	// Reading ne files
	ir = 0;
	iq = 0;
	ih = 0;
	
	for(int i=0; i < nfiles; i++) {
		//strcpy(FILEPATH, ParentDircPath);
		//strcat(FILEPATH, "/");
		//strcat(FILEPATH, nefilename);
		strcpy(FILEPATH, nefilename);
		
		fp[i] = fopen (FILEPATH, "r+");

		//Skips the header rows
		for(j = 0; j < headerrows; j++){
		    fscanf(fp[i], "%*[^\n]\n"); 
		}
	
		// Reads the data into the 2D array line by line
		for(ir = 0; ir < nr_values; ir++){
    		fscanf(fp[i], "%le %le", &r_values[ir], &nedata[iq][ih][ir]);
  		}
  		
  		// Incriments and resets q and h indicies
  		ih++;
  		if (ih >= nH_values) {
  			ih = 0;
  			iq++;
  		} 		
  		if (iq >= nQ_values) iq=0;

	  	fclose(fp[i]);
	  }   
    
    //Saving values into the structure to be retuned
    
    st->Tedata = Tedata;
    st->nedata = nedata;
	st->nr_values = nr_values;
	st->r_values = r_values; 
	st->nQ_values = nQ_values;
	st->Q_values = &Q_values;
	st->nH_values = nH_values;
	st->H_values = &H_values;
	    
	return *st;
}

/* ELECTRON TEMPERATURE FUNCTION: CKC VERSION 1 */
double get_Telec (CKCdata *st, double Q, double rH, double radius) {
    int iQ0, iQ1, ir0, ir1, iH0, iH1; 	// Indicies
    double Qd, Hd, rd; 					// Percent difference
    double T000, T001, T010, T011; 		// T_qhr
    double T100, T101, T110, T111; 		// T_qhr
    double T00, T01, T10, T11; 			// T_hr
    double T0, T1; 						// T_r
    
    int i;
    
    double H = rH;
    int nQ = st->nQ_values;
    int nH = st->nH_values;
    int nr = st->nr_values;
    
    double r = radius * 1e-3; //Converts from m in SUBLIME to km in CKC
       
    double Telec;
 	
    if (Q <= st->Q_values[0]) {
    	iQ1 = 0;
    	iQ0 = 0;
    	Qd = 0;
    	//printf("Case 1: Q = %.2le <= Q%d = %.2le \n", Q, iQ1, st->Q_values[iQ1]);
    }
    else if (Q >= st->Q_values[nQ-1]) {
    	iQ1 = nQ-1;
    	iQ0 = nQ-1;
    	Qd = 0;
    	//printf("Case 2: Q = %.2le >= Q%d = %.2le \n", Q, iQ1, st->Q_values[iQ1]);
    }
    else {
    	for (i = 0; i < nQ; i++) {
			if (Q < st->Q_values[i]) {
				break;
			}
 		}
    	iQ1 = i;
    	iQ0 = iQ1 - 1;
    	Qd = (Q - st->Q_values[iQ0])/(st->Q_values[iQ1] - st->Q_values[iQ0]);
    	//printf("Case 3: Q = %.2le lies between Q%d = %.2le and Q%d = %.2le \n", Q, iQ0, st->Q_values[iQ0], iQ1, st->Q_values[iQ1]);
    }
    
    //printf("Qd = %.2le \n", Qd);
    
    if (H <= st->H_values[0]) {
    	iH1 = 0;
    	iH0 = 0;
    	Hd = 0;
    	//printf("Case A: H = %.2le <= H%d = %.2le \n", H, iH1, st->H_values[iH1]);
    }
    else if (H >= st->H_values[nH-1]) {
    	iH1 = nH-1;
    	iH0 = nH-1;
    	Hd = 0;
    	//printf("Case B: H = %.2le >= H%d = %.2le \n", H, iH1, st->H_values[iH1]);
    }
    else {
    	for (i = 0; i < nH; i++) {
			if (H < st->H_values[i]) {
				break;
			}
 		}
    	iH1 = i;
    	iH0 = iH1 - 1;
    	Hd = (H - st->H_values[iH0])/(st->H_values[iH1] - st->H_values[iH0]);
    	//printf("Case C: H = %.2le lies between H%d = %.2le and H%d = %.2le \n", H, iH0, st->H_values[iH0], iH1, st->H_values[iH1]);
    }
    
    //printf("Hd = %.2le \n", Hd);
        
    if (r <= st->r_values[0]) {
    	ir1 = 0;
    	ir0 = 0;
    	rd = 0;
    	//printf("Case alpha: r = %.2le <= r%d = %.2le \n", r, ir1, st->r_values[ir1]);
    }
    else if (r >= st->r_values[nr-1]) {
    	ir1 = nr-1;
    	ir0 = nr-1;
    	rd = 0;
    	//printf("Case beta: r = %.2le >= r%d = %.2le \n", r, ir1, st->r_values[ir1]);
    }
    else {
    	for (i = 0; i < nr; i++) {
			if (r < st->r_values[i]) {
				break;
			}
 		}
    	ir1 = i;
    	//ir1 = upper_index(r, st->r_values, nr);
    	ir0 = ir1 - 1;
    	rd = (r - st->r_values[ir0])/(st->r_values[ir1] - st->r_values[ir0]);
    	//printf("Case gamma: r = %.2le lies between r%d = %.2le and r%d = %.2le \n", r, ir0, st->r_values[ir0], ir1, st->r_values[ir1]);
    }
    
    //printf("rd = %.2le \n", rd);

    // Identify values to interpolate
    T000 = st->Tedata[iQ0][iH0][ir0];
    T001 = st->Tedata[iQ0][iH0][ir1];
    T010 = st->Tedata[iQ0][iH1][ir0];
    T011 = st->Tedata[iQ0][iH1][ir1];
    T100 = st->Tedata[iQ1][iH0][ir0];
    T101 = st->Tedata[iQ1][iH0][ir1];
    T110 = st->Tedata[iQ1][iH1][ir0];
    T111 = st->Tedata[iQ1][iH1][ir1];
    
    //printf("T000 = %.2le \n", T000);
    //printf("T001 = %.2le \n", T001);
    //printf("T010 = %.2le \n", T010);
    //printf("T011 = %.2le \n", T011);
    //printf("T100 = %.2le \n", T100);
    //printf("T101 = %.2le \n", T101);
    //printf("T110 = %.2le \n", T110);
    //printf("T111 = %.2le \n", T111);
    
    //Begin interpolation
    
    //Thr = T0hr * (1 - Qd) + T1hr * Qd;
    T00 = T000 * (1 - Qd) + T100 * Qd;
    T01 = T001 * (1 - Qd) + T101 * Qd;
    T10 = T010 * (1 - Qd) + T110 * Qd;
    T11 = T011 * (1 - Qd) + T111 * Qd;
    
    //printf("T00 = %.2le \n", T00);
    //printf("T01 = %.2le \n", T01);
    //printf("T10 = %.2le \n", T10);
    //printf("T11 = %.2le \n", T11);
    
    T0 = T00 * (1 - Hd) + T10 * Hd;
    T1 = T01 * (1 - Hd) + T11 * Hd;
    
    //printf("T0 = %.2le \n", T0);
    //printf("T1 = %.2le \n", T1);
    
    Telec = T0 * (1 - rd) + T1 * rd;
    
    /*
    printf("T000 = %.2le \n", T000);
    printf("T100 = %.2le \t T00 = %.2le \n", T100, T00);
    printf("T010 = %.2le \n", T010);
    printf("T110 = %.2le \t T10 = %.2le \t T0 = %.2le \n", T110, T10, T0);
    printf("T001 = %.2le \n", T001);
    printf("T101 = %.2le \t T01 = %.2le \n", T101, T01);
    printf("T011 = %.2le \n", T011);
    printf("T111 = %.2le \t T11 = %.2le \t T1 = %.2le \n", T111, T11, T1);
    */
	return Telec;
}

/* ELECTRON DENSITY FUNCTION: CKC VERSION 1 */
double get_nelec (CKCdata *st, double Q, double rH, double radius) {
    int iQ0, iQ1, ir0, ir1, iH0, iH1; 	// Indicies
    double Qd, Hd, rd; 					// Percent difference
    double n000, n001, n010, n011; 		// n_qhr
    double n100, n101, n110, n111; 		// n_qhr
    double n00, n01, n10, n11; 			// n_hr
    double n0, n1; 						// n_r
        
    int i;
    
    double H = rH;
    int nQ = st->nQ_values;
    int nH = st->nH_values;
    int nr = st->nr_values;
    
    double r = radius * 1e-3; //Converts from m in SUBLIME to km in CKC
    
    double nelec;
 	
    if (Q <= st->Q_values[0]) {
    	iQ1 = 0;
    	iQ0 = 0;
    	Qd = 0;
    	//printf("Case 1: Q = %.2le <= Q%d = %.2le \n", Q, iQ1, st->Q_values[iQ1]);
    }
    else if (Q >= st->Q_values[nQ-1]) {
    	iQ1 = nQ-1;
    	iQ0 = nQ-1;
    	Qd = 0;
    	//printf("Case 2: Q = %.2le >= Q%d = %.2le \n", Q, iQ1, st->Q_values[iQ1]);
    }
    else {
    	for (i = 0; i < nQ; i++) {
			if (Q < st->Q_values[i]) {
				break;
			}
 		}
    	iQ1 = i;
    	iQ0 = iQ1 - 1;
    	Qd = (Q - st->Q_values[iQ0])/(st->Q_values[iQ1] - st->Q_values[iQ0]);
    	//printf("Case 3: Q = %.2le lies between Q%d = %.2le and Q%d = %.2le \n", Q, iQ0, st->Q_values[iQ0], iQ1, st->Q_values[iQ1]);
    }
    
    //printf("Qd = %.2le \n", Qd);
    
    if (H <= st->H_values[0]) {
    	iH1 = 0;
    	iH0 = 0;
    	Hd = 0;
    	//printf("Case A: H = %.2le <= H%d = %.2le \n", H, iH1, st->H_values[iH1]);
    }
    else if (H >= st->H_values[nH-1]) {
    	iH1 = nH-1;
    	iH0 = nH-1;
    	Hd = 0;
    	//printf("Case B: H = %.2le >= H%d = %.2le \n", H, iH1, st->H_values[iH1]);
    }
    else {
    	for (i = 0; i < nH; i++) {
			if (H < st->H_values[i]) {
				break;
			}
 		}
    	iH1 = i;
    	iH0 = iH1 - 1;
    	Hd = (H - st->H_values[iH0])/(st->H_values[iH1] - st->H_values[iH0]);
    	//printf("Case C: H = %.2le lies between H%d = %.2le and H%d = %.2le \n", H, iH0, st->H_values[iH0], iH1, st->H_values[iH1]);
    }
    
    //printf("Hd = %.2le \n", Hd);
        
    if (r <= st->r_values[0]) {
    	ir1 = 0;
    	ir0 = 0;
    	rd = 0;
    	//printf("Case alpha: r = %.2le <= r%d = %.2le \n", r, ir1, st->r_values[ir1]);
    }
    else if (r >= st->r_values[nr-1]) {
    	ir1 = nr-1;
    	ir0 = nr-1;
    	rd = 0;
    	//printf("Case beta: r = %.2le >= r%d = %.2le \n", r, ir1, st->r_values[ir1]);
    }
    else {
    	for (i = 0; i < nr; i++) {
			if (r < st->r_values[i]) {
				break;
			}
 		}
    	ir1 = i;
    	//ir1 = upper_index(r, st->r_values, nr);
    	ir0 = ir1 - 1;
    	rd = (r - st->r_values[ir0])/(st->r_values[ir1] - st->r_values[ir0]);
    	//printf("Case gamma: r = %.2le lies between r%d = %.2le and r%d = %.2le \n", r, ir0, st->r_values[ir0], ir1, st->r_values[ir1]);
    }
    
    //printf("rd = %.2le \n", rd);

    // Identify values to interpolate
    n000 = st->nedata[iQ0][iH0][ir0];
    n001 = st->nedata[iQ0][iH0][ir1];
    n010 = st->nedata[iQ0][iH1][ir0];
    n011 = st->nedata[iQ0][iH1][ir1];
    n100 = st->nedata[iQ1][iH0][ir0];
    n101 = st->nedata[iQ1][iH0][ir1];
    n110 = st->nedata[iQ1][iH1][ir0];
    n111 = st->nedata[iQ1][iH1][ir1];
    
    //printf("n000 = %.2le \n", n000);
    //printf("n001 = %.2le \n", n001);
    //printf("n010 = %.2le \n", n010);
    //printf("n011 = %.2le \n", n011);
    //printf("n100 = %.2le \n", n100);
    //printf("n101 = %.2le \n", n101);
    //printf("n110 = %.2le \n", n110);
    //printf("n111 = %.2le \n", n111);
    
    //Begin interpolation
    n00 = n000 * (1 - Qd) + n100 * Qd;
    n01 = n001 * (1 - Qd) + n101 * Qd;
    n10 = n010 * (1 - Qd) + n110 * Qd;
    n11 = n011 * (1 - Qd) + n111 * Qd;
    //printf("n00 = %.2le \n", n00);
    //printf("n01 = %.2le \n", n01);
    //printf("n10 = %.2le \n", n10);
    //printf("n11 = %.2le \n", n11);
        
    n0 = n00 * (1 - Hd) + n10 * Hd;
    n1 = n01 * (1 - Hd) + n11 * Hd;
    //printf("n0 = %.2le \n", n0);
    //printf("n1 = %.2le \n", n1);
    
    nelec = n0 * (1 - rd) + n1 * rd;
    
    /*
    printf("n000 = %.2le \n", n000);
    printf("n100 = %.2le \t n00 = %.2le \n", n100, n00);
    printf("n010 = %.2le \n", n010);
    printf("n110 = %.2le \t n10 = %.2le \t n0 = %.2le \n", n110, n10, n0);
    printf("n001 = %.2le \n", n001);
    printf("n101 = %.2le \t n01 = %.2le \n", n101, n01);
    printf("n011 = %.2le \n", n011);
    printf("n111 = %.2le \t n11 = %.2le \t n1 = %.2le \n", n111, n11, n1);
    */
	nelec = nelec * 1e6; // Unit conversion: cm^-3 (CKC) -> m^-3 (SUBLIME)

	return nelec;
}