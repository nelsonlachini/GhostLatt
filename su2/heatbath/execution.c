#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <string.h>

#include "ran0.h"

#include "utilities.c"
#include "algebra.c"
#include "thermalization_hb.c"

#define MAXIMUM_NAME_LENGTH 50

//GLOBAL VARIABLES
int N;
double beta;

void main(){
	//LOCAL PARAMETERS	
	double * lattice;	
	long * global_seed = malloc(sizeof(long));	
	int i, j, N_cf , N_g;	
	float initial_order;	

	N = 8;
	N_cf = 5;
	N_g = 15;
	beta = 2.45;
	initial_order = 0;					//0 to all links initially equal to identity	
	
	//EXECUTION
	lattice = malloc(sizeof(double)*N*N*N*N*4*4);	
		
	initl( lattice , initial_order , global_seed);
	thermalizeLattice(lattice , N_cf , global_seed);	

	//GARBAGE
	free(lattice);		
	free(global_seed);
}