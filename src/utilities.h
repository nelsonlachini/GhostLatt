#ifndef UTILITIES_H
#define UTILITIES_H

#include <complex.h>
#include <stdio.h>

//////////////////////////////ran0 implementation

float ran0(long *idum);

////////////////////////////// SU2 MATRIX/VECTORS METHODS

double * getelv(double * vector , int i);

double complex getelm(double * vector , int i, int j);

void printv(double * vector);

void printc(double complex z);

void printr(double r);

void fprintc( FILE * file , double complex z);

void printm(double * A);

void fprintm(FILE * file , double * A);

double * getU(double * lattice, int t , int x , int y, int z , int mi);

double * getUdim(double * lattice, int t , int x , int y, int z , int mi, int N);

void printLattice(double * lattice);

void copyv(double * out , double * in);

double * randSU2v(double * vector ,long * seed , double parameter);


////////////////////////////// GENERAL MATRIX FUNCTIONS WITH DOUBLE VARIABLES

double * getelmr(double * matrix , int i , int j , int dim);

void printmr(double * A , int dim);

void copymr(double * Mout , double * Min , int dim);

void copyvr(double * vout, double * vin , int dim);

double * getelvr(double * vr , int a , int t , int x , int y , int z);

void printvr(double * V , int dim);

////////////////////////////// GENERAL MATRIX FUNCTIONS WITH DOUBLE COMPLEX VARIABLES

double complex * getelmc(double complex * M , int i , int j , int dim);

void printvc(double complex * V , int dim);

void copyvc(double complex * Vout , double complex * Vin, int dim);

void copymc(double complex * Mout , double complex * Min , int dim);

void printmc(double complex * A , int dim);

double complex * getelvc(double complex * vc , int a , int t , int x , int y , int z);

void converttov(double * vout, double complex * vin);

////////////////////////////// OTHER STUFF

void setquadv(int * amu, int mu);

void initl(double * lattice , float parameter);

void reunitl(double * lattice);

void reunitv(double * vector);

double modulus(double a, double b);

double getStep( double coord , double dcoord );

double getStepT( double coord , double dcoord );

void copyl(double * Lout , double * Lin);

double * getg(double * g, int t , int x , int y, int z);

void printvcbycolor(double complex * v, int a);

void tseed();

void progress_panel(int i, int total);

void save_lattice(double * lattice, char * file_name);

void save_vector(double complex * vector, char * file_name);

void save_rvector(double * vector, char * file_name);

void getName(char * out, char * header , double num);

void load_lattice(double * lattice , char * file_name);

void load_rvector(double * vector , char * file_name);

int compareLink(double * link1, double * link2);

int compareLattice(double * lattice1, double * lattice2);

int comparevc(double complex * v1, double complex * v2, int dim);

int comparevr(double * v1, double * v2, int dim);

void reescaleGaugeField(double * out, double * in, double scale);

void reescaleLinks(double * out, double * in, double scale);

double compareVectorToPW(double * eigenvector, double pPW);

double zeroMomentumTransform(double * vector);

int isLatticeUnitary(double * lattice);

int isLinkUnitary(double * link);

int isLatticeNaN(double * lattice);

double max(double * array, int dim);

double min(double * array, int dim);

double sortDouble(double i , double j);

void printpos(double t, double x, double y , double z, double mi);

#endif
