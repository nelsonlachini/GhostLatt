#ifndef THERMALIZATION_HB_H
#define THERMALIZATION_HB_H

#include "algebra.h"

double updateStaple(double * lattice , int t, int x, int y, int z, int mi , double * _staple);

void su2Bath( double * X , double a );

void microCanonicalStep(double * link, double * _staple);

void microCanonicalUpdate(double * lattice);

void heatBathStep(double * link, double a, double * Vdagger);

void heatbathUpdate(double * lattice);

void overrelaxationStepSU2(double * lattice , int t, int x, int y, int z, int mi, double * Vdagger);

double updateLattice(double * lattice, int N_hb, int N_mc);

double thermalizeLattice(double * lattice , int n , int n_hb , int n_over);

#endif
