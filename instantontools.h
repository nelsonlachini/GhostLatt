#ifndef INSTANTON_TOOLS_H
#define INSTANTON_TOOLS_H

#include <stdio.h>

//////////////////////////tools

void getImagv(double * Mout, double * vin);

void plaquetteProd(double * Vout, double * lattice , int t, int x, int y, int z, int mu, int nu);

void plaquetteProdG(double * Vout, double * lattice , int t, int x, int y, int z, int mu, int nu);

void cloverProd(double * Vout, double * lattice , int t, int x, int y, int z, int mu, int nu);

double integrate(double * indens);

void integrate1d(double * outdenst, double * indens);

void integrate2d(double * outdenst, double * indens);

///////////////////////////smoothings
void coolingStep(double * lattice);

double coolLattice(double * lattice, int Ncool);

double coolLatticeBorders(double * lattice, int Ncool);

void APEsmearStep(double * lattice, double alphaape);

double APEsmearLattice(double * lattice, double alphaape, int Ncool);

/////////////////////action definitions
double calcSwilson(double * lattice, double * sdens);

/////////////////////topological charge field-theoretical definitions
double calcQnaive(double * lattice, double * qdens);

double calcQnaiveimag(double * lattice , double * qdens);

//REPAIR DUE TO P.B.C TO BE MADE
double calcQcloverimag(double * lattice , double * qdens);

double calcQnaivesym(double * lattice,double * qdens);

//AUXILIAR

void saveChargeProfiles(double * Qdensity, char file_name[50], FILE * f);

double continuum1dChargeDensity(double rho, double x, double x0);

#endif
