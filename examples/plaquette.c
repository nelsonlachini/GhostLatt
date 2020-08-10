#include <stdlib.h>
#include <complex.h>
#include <stdio.h>

//maximum set of includes
#include <src/global.h>
#include "src/algebra.h"
#include "src/utilities.h"
#include <src/thermalization_hb.h>
#include <src/measurement.h>

// global variables (to do: change that to local)
const int N = 8;

const int Nt = N;
const int totalV = N*N*N*Nt;
const int spatialV = N*N*N;
const int dimLattice = totalV*4*4;
const int colorV = totalV*3;
const int N2 = N*N;

long * global_seed;

int main(){
    tseed();			//setting global seed

    //LOCAL PARAMETERS
    double beta = 2.2;
    double initial_order = 0;		// 0 is cold start
    int N_therm = 1000;		        // #thermalization sweeps
    int N_hb 	= 1;		        // #heat-bath steps per sweep
    int N_mic 	= N/2;	            // #overrelaxation steps per sweep

    LatticeLinkSU2 * lattice = newLatticeLinkSU2(N,Nt); 

    initl(lattice , initial_order);			//initialize the lattice

    // HOR thermalization and plaquette  measurement
    for(unsigned int i=0 ; i<N_therm ; i++){
    	thermalizeLattice(lattice, beta, 1, N_hb, N_mic);
	    printf("\rHOR Sweep %d : plaquette = %lf",i,wilsonLoopMeasureSym(lattice,1,1));	
    }
    
    free(lattice->U);

    return(0);
}
