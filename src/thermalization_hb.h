#ifndef THERMALIZATION_HB_H
#define THERMALIZATION_HB_H

#include "algebra.h"

void updateStaple(LatticeSU2 * lattice , int t, int x, int y, int z, int mi , double * _staple);

void su2Bath( double * X , double a , double beta);

void microCanonicalStep(double * link, double * _staple);

void microCanonicalUpdate(LatticeSU2 * lattice);

void heatBathStep(double * link, double a, double * Vdagger , double beta);

void heatbathUpdate(LatticeSU2 * lattice , double beta);

void overrelaxationStepSU2(LatticeSU2 * lattice , int t, int x, int y, int z, int mi, double * Vdagger);

double updateLattice(LatticeSU2 * lattice , double beta , int N_hb, int N_mc);

double thermalizeLattice(LatticeSU2 * lattice , double beta , int n , int n_hb , int n_over);

#endif
