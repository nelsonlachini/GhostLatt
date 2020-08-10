#include <stdlib.h>
#include <complex.h>
#include <stdio.h>

#include <src/global.h>
#include "src/algebra.h"
#include "src/utilities.h"
#include <src/thermalization_hb.h>
#include <src/measurement.h>
#include <src/gaugefix.h>
#include <src/statistics.h>
#include <src/inverters.h>

// global variables (to do: change that to local)
const int N = 8;

const int Nt = N;
const int totalV = N*N*N*Nt;
const int spatialV = N*N*N;
const int dimLattice = totalV*4*4;
const int colorV = totalV*3;
const int N2 = N*N;

FILE * f;
long * global_seed;
double THERM_TIME=0, FIX_TIME=0, CG_TIME=0;

int main(){
    clock_t dtime = clock();
    tseed();

	//LOCAL PARAMETERS
    int i,j;
    double beta = 2.7;
	int N_therm = 2200;
	int N_cor = 200;
	int N_cf = 400;
    int N_hb = 1;                                   //hb steps per N_cor
    int N_mic = N/2;                                //overrelaxation steps per N_cor
    int N_bootstrap = 3000;
	double initial_order = 0e0;		        	    //0 to all links initially equal to identity
    char file_name[MAXIMUM_NAME_LENGTH], aux[15];
    
    LatticeLinkSU2 * lattice = newLatticeLinkSU2(N,Nt);                     //lattice alocation
    LatticeColorVectorReal * eigenvector_guess = newLatticeColorVectorReal(lattice->N,lattice->Nt);

	double * lambda1 = malloc(sizeof(double)*N_cf);

    //uncorrelation and measurement steps
    LatticeGaugeSU2 * g = newLatticeGaugeSU2(lattice->N,lattice->Nt);   //allocating gauge transformation
    setidLatticeGaugeSU2(g);                                       //initializing gauge transformation

    double e2tol = 1e-14, r_CG=0, eigen_tol=1e-6 , cg_tol=1e-10;                           // gauge fixing precision
    int iCG;
    double p_stoch = 0.8; //= calibrate_stoc(lattice , 1e-10);    // choosing or calibrating gauge fixing parameter
    
    double pmin[4] = {0e0,0e0,0e0,1e0/lattice->N};

    initl(lattice , initial_order);                                              //cold initialization
    THERM_TIME += thermalizeLattice(lattice, beta, N_therm, N_hb, N_mic);        //thermalization sweeps

    //measurements on the go
    for(i=0;i<N_cf;i++){
        printf("\n########## Measurement step %d\n",i);
        THERM_TIME += thermalizeLattice(lattice, beta, N_cor, N_hb, N_mic);      //N_cor uncorrelation steps of HOR type

        FIX_TIME += fixLatticeStoch(lattice, lattice, g, p_stoch, e2tol);

        printf("\nPower Method execution:");
        setplanewaveallcolorsLatticeColorVectorReal(eigenvector_guess, pmin);
        // compute smallest eigenvalue and eigenvector (output in 2nd and 3rd arguments)
		CG_TIME += rsmallest_eigen_cg(lattice, &lambda1[i], eigenvector_guess, eigen_tol , eigenvector_guess, &r_CG , &iCG, cg_tol);
    }

    double lambda_Lap = 4*pow(sin(M_PI/lattice->N),2);
    //bootstrap analysis
    double lambda1Boot[2];
    bootstrapSimple(&lambda1Boot[0], &lambda1[0], N_cf, N_cf, N_bootstrap);
    printf("\nlambda_1=%lf+-%lf" , lambda1Boot[0] , lambda1Boot[1]);
    printf("\nlambda_1/lambda_Lap=%lf+-%lf" , lambda1Boot[0]/lambda_Lap , lambda1Boot[1]/lambda_Lap);

    ////////////////////////////////PRINTING DATA
    sprintf(file_name,"%s","smallestEigenFP");
    sprintf(aux, "_N%d", N);
    strcat(file_name, aux);
    sprintf(aux, "_B%.2f", beta);
    strcat(file_name, aux);
    printf("\n\n#### Saving results to %s",file_name);
    /////////////////////////////////////////////
    f = fopen(file_name, "w");
	fprintf(f, "#N = %d\n", N);
	fprintf(f, "#BETA = %.2lf\n", beta);
	fprintf(f, "#N_cf = %d\n", N_cf);
	fprintf(f, "#N_therm= %d\n", N_therm);
	fprintf(f, "#N_cor = %d\n", N_cor);
    fprintf(f, "#N_hb = %d\n", N_hb);
    fprintf(f, "#N_mic = %d\n", N_mic);
	fprintf(f, "#INITIAL ORDER OF LINKS = %.2lf\n" , initial_order);
    fprintf(f, "#THERM_TIME: %.2lf hrs\n", THERM_TIME/3600.0);
    fprintf(f, "#FIX_TIME: %.2lf hrs\n", FIX_TIME/3600.0);
    fprintf(f, "#CG_TIME: %.2lf hrs\n", CG_TIME/3600.0);
	fprintf(f, "#TOTAL EXECUTION TIME: %.2lf hrs\n", ((double)(clock() - dtime))/(CLOCKS_PER_SEC*3600.0));
    fprintf(f,"#####################PROPAGATOR SECTION#######################\n");
    fprintf(f,"lambda_1=%lf+-%lf\n" , lambda1Boot[0] , lambda1Boot[1]);
    fprintf(f,"lambda_1/lambda_Lap=%lf+-%lf\n" , lambda1Boot[0]/lambda_Lap , lambda1Boot[1]/lambda_Lap);
    fprintf(f,"#####################SAMPLES DUMP#######################\n");
    fprintf(f,"lambda1");
    for(j=0;j<N_cf;j++){
        fprintf(f,"%lf\n",lambda1[j]);
    }
    fclose(f);

    freeLatticeLinkSU2(lattice);
    freeLatticeGaugeSU2(g);
    freeLatticeColorVectorReal(eigenvector_guess);
    free(lambda1);

    return(0);
}