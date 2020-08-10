#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "algebra.h"
#include "global.h"

double wilsonLoopMeasureSym(LatticeLinkSU2 * lattice , int X, int Y);

double wilsonLoopMeasureAsym(LatticeLinkSU2 * lattice , int R, int T);

double polyakovLoopMeasure(LatticeLinkSU2 * lattice);

double gluonP(LatticeLinkSU2 * lattice, double kz);

double measure_ghostp(LatticeLinkSU2 * lattice, double *** G, double * kz, int kSize, int j_Ncf);

double measure_gluonp(LatticeLinkSU2 * lattice, double *** D, double * kz, int kSize, int j_Ncf);

void updateStapleMeasurement(LatticeLinkSU2 * lattice , int t, int x, int y, int z, int mi , double * _staple);

void spatialSmearStep(LatticeLinkSU2 * outlattice, LatticeLinkSU2 * inlattice, double alpha_smear);

void spatialSmearing3(LatticeLinkSU2 * outlattice, LatticeLinkSU2 * inlattice, int n_smear, double alpha_smear);

#endif
