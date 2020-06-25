#ifndef GAUGE_FIX_H
#define GAUGE_FIX_H

#include "global.h"

//////////////////////////auxiliars
void initg(double * g);

double applyFix(LatticeSU2 * Uout, LatticeSU2 * Uin , double * g);

void copyg(double * gout, double * gin);

double su2distance(LatticeSU2 * U, LatticeSU2 * W);

//SERIAL VERSIONS

double los_alamos(LatticeSU2 * U , double * g , double * e2);

double cornell(LatticeSU2 * U , double * g, double * e2 , double alpha);

double overrelaxation(LatticeSU2 * U , double * g, double * e2 , double omega);

double stochastic(LatticeSU2 * U , double * g , double * e2, double p);

// double fourier(double * U , double * g , double * e2, double alpha);

//AUTOMATIZATION
double fixLatticeStoch(LatticeSU2 * lout, LatticeSU2 * lin, double * g, double p_stoch, double e2tol);

double fixLatticeLosalamos(LatticeSU2 * lout, LatticeSU2 * lin, double * g , double e2tol);

double fixLatticeOver(LatticeSU2 * lout, LatticeSU2 * lin, double * g, double omega, double e2tol);

//AUTO CALIBRATION

double calce2(LatticeSU2 * U);

double calibrate_over(LatticeSU2 * lattice , double tol);

double calibrate_cornell(LatticeSU2 * lattice , double tol, double step_size, double lower_alpha, double upper_alpha);

double calibrate_stoc(LatticeSU2 * lattice , double tol);

// double * getQ(int ni , int xni);

// double * getAvgQ(int ni);

// void calcAvgQ(int ni , double * sum);

// void calcQni( int ni  , double * lattice);

//double calce6(double * lattice , double * g);

//####################################################IN TEST
void reunitg(double * g);

double calcEps(LatticeSU2 * lattice , double * g);

double calcFunctional(double * U);

double * getLatticeReal(double * real , int t , int x , int y , int z);

double latticeAvgReal(double * var);

#endif
