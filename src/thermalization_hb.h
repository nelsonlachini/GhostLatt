#ifndef THERMALIZATION_HB_H
#define THERMALIZATION_HB_H

#include "algebra.h"

void updateStaple(LatticeLinkSU2 * lattice , int t, int x, int y, int z, int mi , double * _staple);

void su2Bath( double * X , double a , double beta);

void microCanonicalStep(double * link, double * _staple);

void microCanonicalUpdate(LatticeLinkSU2 * lattice);

void heatBathStep(double * link, double a, double * Vdagger , double beta);

void heatbathUpdate(LatticeLinkSU2 * lattice , double beta);

void overrelaxationStepSU2(LatticeLinkSU2 * lattice , int t, int x, int y, int z, int mi, double * Vdagger);

double updateLattice(LatticeLinkSU2 * lattice , double beta , int N_hb, int N_mc);

double thermalizeLattice(LatticeLinkSU2 * lattice , double beta , int n , int n_hb , int n_over);

#endif
