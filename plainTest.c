
//Complete compile command: 
//gcc-8 plainTest.c utilities.c algebra.c thermalization_hb.c measurement.c 
//  inverters.c gaugefix.c instantontools.c nrutil.c statistics.c hoek_custom.c
//  fox.c center.c -o temp -lm

#include <stdlib.h>
#include <complex.h>
#include <stdio.h>

#include "global.h"
#include "thermalization_hb.h"
#include "measurement.h"
#include "gaugefix.h"
#include "inverters.h"
#include "instantontools.h"
#include "nrutil.h"
#include "statistics.h"
#include "hoek_custom.h"
#include "fox.h"
#include "center.h"

// global variables (need to change that)
const int N = 8;
const double beta = 2.45;

const int Nt = N;
const int totalV = N*N*N*Nt;
const int spatialV = N*N*N;
const int dimLattice = totalV*4*4;
const int colorV = totalV*3;
const int N2 = N*N;

long * global_seed;

int main(){
    tseed();

	//LOCAL PARAMETERS
	int N_therm = 10;
    int N_hb = 1;                               //hb steps per N_cor
    int N_mic = N/2;                            //overrelaxation steps per N_cor
    
	double initial_order = 0e0;					//0 to all links initially equal to identity

    double * thermlattice = malloc(sizeof(double)*dimLattice);

    initl(thermlattice , initial_order);
    thermalizeLattice(thermlattice, N_therm, N_hb, N_mic);
    
    free(thermlattice);

    return 1;
}
