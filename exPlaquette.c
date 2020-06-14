
#include <stdlib.h>
#include <complex.h>
#include <stdio.h>

//maximum set of includes
#include "src/global.h"
#include "src/thermalization_hb.h"
#include "src/measurement.h"
#include "src/gaugefix.h"
#include "src/inverters.h"
#include "src/instantontools.h"
#include "src/nrutil.h"
#include "src/statistics.h"
#include "src/hoek_custom.h"
#include "src/fox.h"
#include "src/center.h"

// global variables (to do: change that to local)
const int N = 8;
const double beta = 2.2;

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
    int N_therm = 1000;		// #thermalization sweeps
    int N_hb 	= 1;		// #heat-bath steps per sweep
    int N_mic 	= N/2;	        // #overrelaxation steps per sweep
    
    double initial_order = 0;			//0 is cold start

    double * lattice = malloc(sizeof(double)*dimLattice);	//allocating lattice (to do: encapsulate it)

    initl(lattice , initial_order);			//initialize the lattice

    // HOR thermalization while measuring plaquette
    for(unsigned int i=0 ; i<N_therm ; i++){
    	thermalizeLattice(lattice, 1, N_hb, N_mic);
	    printf("\rHOR Sweep %d : plaquette = %lf",i,wilsonLoopMeasureSym(lattice,1,1));	
    }
    
    free(lattice);

    return 0;
}
