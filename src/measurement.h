#ifndef MEASUREMENT_H
#define MEASUREMENT_H

double wilsonLoopMeasureSym(double * lattice , int X, int Y);

double wilsonLoopMeasureAsym(double * lattice , int R, int T);

double polyakovLoopMeasure(double * lattice);

double measure_ghostp(double * lattice, double *** G, double * kz, int kSize, int j_Ncf);

double measure_gluonp(double * lattice, double *** D, double * kz, int kSize, int j_Ncf);

void spatialSmearing3(double * outlattice, double * inlattice, int n_smear, double alpha_smear);

#endif
