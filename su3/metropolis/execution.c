#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "mcsimulation.c"

////////////////////////////declaration of global variables
int N_cf;
int N;
int N_hit;
double beta;
double eps;
double complex * pool;

void mcAverage(long * seed ){
	int i;	
	printf("\nW(1,1)");
	for(i = 0 ; i < N_cf ; i++){
		if (i % 200 == 0){
			createPoolSU3(seed);
		}
		updateLattice(seed);				
		printf("\n%lf", wilsonLoopMeasure(1,1));
	}
	printf("\n\nBETA = %lf",beta);
}

void main(){
	long * global_seed;	
	float initial_order;
	
	N = 8;
	beta = 5.5;
	N_cf = 10;
	N_hit = 10;
	eps = 0.24;
	initial_order = 0;					//0 to all links initially equal to identity

	global_seed = genSeed();
	
	lattice = malloc(sizeof(double complex)*N*N*N*N*SP_DIM*D*D);
	pool = malloc(sizeof(double complex)*2*POOL_SIZE*D*D);
	
	initLattice(0);
	createPoolSU3(global_seed);
	printMatrix(getFromPool(0));
	printMatrix(getFromPool(1));
	printMatrix(getFromPool(2));
	printMatrix(getFromPool(3));
	
	//mcAverage(global_seed);

	free(global_seed);
	free(lattice);	
}