#ifndef GAUGE_FIX_H
#define GAUGE_FIX_H

#include "global.h"

//////////////////////////auxiliars
void setidLatticeGaugeSU2(LatticeGaugeSU2 * g);

void reunitLatticeGaugeSU2(LatticeGaugeSU2 * g);

double calcEps(LatticeLinkSU2 * lattice , LatticeGaugeSU2 * g);

double applyFix(LatticeLinkSU2 * Uout, LatticeLinkSU2 * Uin , LatticeGaugeSU2 * g);

void copyg(LatticeGaugeSU2 * gout, LatticeGaugeSU2 * gin);

double su2distance(LatticeLinkSU2 * U, LatticeLinkSU2 * W);

//SERIAL VERSIONS

double los_alamos(LatticeLinkSU2 * U , LatticeGaugeSU2 * g , double * e2);

double cornell(LatticeLinkSU2 * U , LatticeGaugeSU2 * g, double * e2 , double alpha);

double overrelaxation(LatticeLinkSU2 * U , LatticeGaugeSU2 * g, double * e2 , double omega);

double stochastic(LatticeLinkSU2 * U , LatticeGaugeSU2 * g , double * e2, double p);

// double fourier(double * U , double * g , double * e2, double alpha);

//AUTOMATIZATION
double fixLatticeStoch(LatticeLinkSU2 * lout, LatticeLinkSU2 * lin, LatticeGaugeSU2 * g, double p_stoch, double e2tol);

double fixLatticeLosalamos(LatticeLinkSU2 * lout, LatticeLinkSU2 * lin, LatticeGaugeSU2 * g , double e2tol);

double fixLatticeOver(LatticeLinkSU2 * lout, LatticeLinkSU2 * lin, LatticeGaugeSU2 * g, double omega, double e2tol);

//AUTO CALIBRATION

double calce2(LatticeLinkSU2 * U);

double calibrate_over(LatticeLinkSU2 * lattice , double tol);

double calibrate_cornell(LatticeLinkSU2 * lattice , double tol, double step_size, double lower_alpha, double upper_alpha);

double calibrate_stoc(LatticeLinkSU2 * lattice , double tol);

//####################################################IN TEST

double calcFunctional(LatticeLinkSU2 * lattice);

#endif
