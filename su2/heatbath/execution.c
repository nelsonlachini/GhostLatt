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

//GLOBAL VARIABLES
int N;
double beta;

void main(){
	//LOCAL PARAMETERS
	double * lattice;
	long * global_seed = malloc(sizeof(long));
	int i, j, N_cf , N_g;
	float initial_order;							//this parameter varies from 0 to 1; 0 correspond to all links 
													//	initialized as identity; 1 correspond to the other extreme

	N = 8;
	N_cf = 20;
	N_g = 15;
	beta = 2.45;
	initial_order = 0;								//0 to all links initially equal to identity
	
	//EXECUTION
	lattice = malloc(sizeof(double)*N*N*N*N*4*4);
		
	initl( lattice , initial_order , global_seed);
	thermalizeLattice(lattice , N_cf , global_seed);
	
	//GARBAGE
	free(lattice);
	free(global_seed);
}