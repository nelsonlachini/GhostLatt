#ifndef UTILITIES_H
#define UTILITIES_H

#include <complex.h>
#include <stdio.h>

#include "global.h"

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

double * getLink(LatticeLinkSU2 * lattice, int t , int x , int y, int z , int mi);

double * getLinkdim(LatticeLinkSU2 * lattice, int t , int x , int y, int z , int mi, int N);

void printLattice(LatticeLinkSU2 * lattice);

void copyv(double * out , double * in);

void randSU2v(double * vector ,long * seed , double parameter);

////////////////////////////// GENERAL MATRIX FUNCTIONS WITH DOUBLE VARIABLES

double * getelmr(double * matrix , int i , int j , int dim);

void printmr(double * A , int dim);

void copymr(double * Mout , double * Min , int dim);

void copyvr(double * vout, double * vin, int dim);

void copyLatticeColorVectorReal(LatticeColorVectorReal * vout, LatticeColorVectorReal * vin);

void copyLatticeColorVectorComplex(LatticeColorVectorComplex * vout, LatticeColorVectorComplex * vin);

double * getelvr(double * vr , int a , int t , int x , int y , int z);

double * getLatticeColorVectorReal(LatticeColorVectorReal * vr , int a , int t , int x , int y , int z);

void printvr(double * V , int dim);

////////////////////////////// GENERAL MATRIX FUNCTIONS WITH DOUBLE COMPLEX VARIABLES

double complex * getelmc(double complex * M , int i , int j , int dim);

void printvc(double complex * V , int dim);

void copyvc(double complex * Vout , double complex * Vin, int dim);

void copymc(double complex * Mout , double complex * Min , int dim);

void printmc(double complex * A , int dim);

double complex * getelvc(double complex * vc , int a , int t , int x , int y , int z);

double complex * getLatticeColorVectorComplex(LatticeColorVectorComplex * vc , int a , int t , int x , int y , int z);

void converttov(double * vout, double complex * vin);

////////////////////////////// OTHER STUFF

void setUnitVector(int * amu, int mu);

void initl(LatticeLinkSU2 * lattice , float parameter);

void reunitLatticeLinkSU2(LatticeLinkSU2 * lattice);

void reunitv(double * vector);

double modulus(double a, double b);

double getStep( double coord , double dcoord );

double getStepT( double coord , double dcoord );

void copyLatticeLinkSU2(LatticeLinkSU2 * Lout , LatticeLinkSU2 * Lin);

double * getg(LatticeGaugeSU2 * g, int t , int x , int y, int z);

void printvcbycolor(double complex * v, int a);

void tseed();

void progress_panel(int i, int total);

void saveLatticeLinkSU2(LatticeLinkSU2 * lattice, char * file_name);

void saveLatticeColorVectorComplex(LatticeColorVectorComplex * vector, char * file_name);

void saveLatticeColorVectorReal(LatticeColorVectorReal * vector, char * file_name);

void getName(char * out, char * header , double num);

void loadLatticeLinkSU2(LatticeLinkSU2 * lattice , char * file_name);

void loadLatticeColorVectorReal(LatticeColorVectorReal * vector , char * file_name);

void loadLatticeColorVectorComplex(LatticeColorVectorComplex * vector , char * file_name);

int compareLink(double * link1, double * link2);

int compareLattice(LatticeLinkSU2 * lattice1, LatticeLinkSU2 * lattice2);

int compareLatticeColorVectorComplex(LatticeColorVectorComplex * v1, LatticeColorVectorComplex * v2);

int compareLatticeColorVectorReal(LatticeColorVectorReal * v1, LatticeColorVectorReal * v2);

void reescaleGaugeField(LatticeLinkSU2 * out, LatticeLinkSU2 * in, double scale);

void reescaleLinks(LatticeLinkSU2 * out, LatticeLinkSU2 * in, double scale);

double compareVectorToPW(LatticeColorVectorReal * eigenvector, double pPW);

double zeroMomentumTransform(double * vector);

int isLatticeUnitary(LatticeLinkSU2 * lattice);

int isLinkUnitary(double * link);

int isLatticeNaN(LatticeLinkSU2 * lattice);

double max(double * array, int dim);

double min(double * array, int dim);

double sortDouble(double i , double j);

void printpos(double t, double x, double y , double z, double mi);

LatticeLinkSU2* newLatticeLinkSU2(int _N, int _Nt);
LatticeColorVectorReal* newLatticeColorVectorReal(int _N, int _Nt);
LatticeColorVectorComplex* newLatticeColorVectorComplex(int _N, int _Nt);
LatticeGaugeSU2* newLatticeGaugeSU2(int _N, int _Nt);

void freeLatticeLinkSU2(LatticeLinkSU2 * l);
void freeLatticeColorVectorReal(LatticeColorVectorReal * v);
void freeLatticeColorVectorComplex(LatticeColorVectorComplex * v);
void freeLatticeGaugeSU2(LatticeGaugeSU2 * l);

#endif
