#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mcsimulation_m.c"

#define MAXIMUM_NAME_LENGTH 50

////////////////////////////declaration of global variables
int N_cf;
int N;
int N_hit;
double beta;
double eps;

void mcAverage(long * seed ){	
	int i;
	printf("\n W(1,1)");
	for(i = 0 ; i < N_cf ; i++){
		updateLattice(seed);							
		printf("\n %lf",wilsonLoopMeasure(1,1));		
	}
	printf("\n\n BETA = %lf \n" , beta);
}

void main(){
	long * global_seed;	
	int j;	
	float initial_order;	
		
	N=2;
	beta=2.3;
	eps = 0.2;
	N_cf = 10;
	N_hit = 10;
	initial_order = 1.;
		
	lattice = malloc(sizeof(double complex)*N*N*N*N*SP_DIM*D*D);
	pool = malloc(sizeof(double complex)*2*POOL_SIZE*D*D);
	
	global_seed = genSeed();
	
	createPool(global_seed);
	initLattice( initial_order , global_seed);
	mcAverage(global_seed);	
			
	free(global_seed);
	free(lattice);	
}