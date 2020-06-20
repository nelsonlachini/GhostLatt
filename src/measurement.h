#ifndef MEASUREMENT_H
#define MEASUREMENT_H

#include "algebra.h"
#include "global.h"

double wilsonLoopMeasureSym(LatticeSU2 * lattice , int X, int Y);

double wilsonLoopMeasureAsym(LatticeSU2 * lattice , int R, int T);

double polyakovLoopMeasure(LatticeSU2 * lattice);

double gluonP(LatticeSU2 * lattice, double kz);

double measure_ghostp(LatticeSU2 * lattice, double *** G, double * kz, int kSize, int j_Ncf);

double measure_gluonp(LatticeSU2 * lattice, double *** D, double * kz, int kSize, int j_Ncf);

void updateStapleMeasurement(LatticeSU2 * lattice , int t, int x, int y, int z, int mi , double * _staple);

void spatialSmearStep(LatticeSU2 * outlattice, LatticeSU2 * inlattice, double alpha_smear);

void spatialSmearing3(LatticeSU2 * outlattice, LatticeSU2 * inlattice, int n_smear, double alpha_smear);

#endif
