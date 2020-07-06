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

// global variables (to do: change that to local)
const int N = 8;

const int Nt = N;
const int totalV = N*N*N*Nt;
const int spatialV = N*N*N;
const int dimLattice = totalV*4*4;
const int colorV = totalV*3;
const int N2 = N*N;

long * global_seed;
FILE * f;

int main(){
    clock_t dtime = clock();
    double therm_time;
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
    
    LatticeSU2 * lattice = newLatticeSU2(N,Nt);                     //lattice alocation

    initl(lattice , initial_order);                                //cold initialization
    thermalizeLattice(lattice, beta, N_therm, N_hb, N_mic);        //thermalization sweeps

    int kListSize = (int)(0.5*N); //half the lattice
    double * kList = malloc(sizeof(double)*kListSize);
    double ** ghostProp;
    ghostProp = (double **)malloc(sizeof(double *)*kListSize);
    for(i=0;i<kListSize;i++){
        ghostProp[i] = (double *)malloc(sizeof(double)*N_cf);
    }

    //uncorrelation and measurement steps
    double * g = malloc(sizeof(double)*totalV*4);   //allocating gauge transformation
    initg(g);                                       //initializing gauge transformation
    double e2tol = 1e-14;                           // gauge fixing precision
    double p_stoch = 0.8; //= calibrate_stoc(lattice , 1e-10);    // choosing or calibrating gauge fixing parameter

    for(i=0;i<N_cf;i++){
        printf("\n########## Measurement step %d\n",i);
        thermalizeLattice(lattice, beta, N_cor, N_hb, N_mic);    //N_cor uncorrelation steps of HOR type

        fixLatticeStoch(lattice, lattice, g, p_stoch, e2tol);

        printf("\nMeasurements:");
        measure_ghostp(lattice, &ghostProp , kList, kListSize, i);
    }

    //bootstrap analysis
    double ** ghostPropBoot = malloc(sizeof(double *)*kListSize);
    for(i=0;i<kListSize;i++){
        ghostPropBoot[i] = malloc(sizeof(double)*2);
        bootstrapSimple(ghostPropBoot[i], ghostProp[i], N_cf, N_cf, N_bootstrap);
        printf("\nG(%.3lf)=%lf+-%lf" , kList[i] , ghostPropBoot[i][0] , ghostPropBoot[i][1]);
    }

    ////////////////////////////////PRINTING DATA
    sprintf(file_name,"%s","ghostPropagator");
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
	fprintf(f, "#TOTAL EXECUTION TIME: %.2lf hrs\n", ((double)(clock() - dtime))/(CLOCKS_PER_SEC*3600.0));
    fprintf(f,"#####################PROPAGATOR SECTION#######################\n");
    fprintf(f,"#k/N , Lattice momentum , D(k) , stdD \n");
    for(i=0;i<kListSize;i++){
        fprintf(f,"%.3lf %.3f %lf %lf \n", kList[i] , 2*sin(M_PI*kList[i]) , ghostPropBoot[i][0] , ghostPropBoot[i][1]);
    }
    fprintf(f,"#####################SAMPLES DUMP#######################\n");
    fprintf(f,"\n");
    for(i=0;i<kListSize;i++){
        fprintf(f,"D(%.3lf)|",kList[i]);
        for(j=0;j<N_cf;j++){
            fprintf(f,"%lf|",ghostProp[i][j]);
        }
        fprintf(f,"\n");
    }
    fclose(f);

    free(g);
    free(lattice->U);

    return(0);
}
