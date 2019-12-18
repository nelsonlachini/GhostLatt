#ifndef CENTER_H
#define CENTER_H

double largest_eigen(double * matrix, int dim,  double * lambda1, double  * eigen_out);

void constructD(double * lattice, double * D, int t, int x, int y, int z);

double directMCGOverrelaxationSweep(double * lattice, double * g, double omega);

double directMCGStochasticSweep(double * lattice, double * g, double p_stoch);

double calcR(double * lattice);

void centerProjection(double * vortexlattice, double * MCGlattice);

void centerRemoval(double * MCGlattice, double * vortexlattice);

double MCGOver(double * outlattice, double * inlattice, /* double * g, */ double omega, double Rtol);

double MCGStoch(double * outlattice, double * inlattice, /* double * g, */ double p_stoch, double Rtol);

void centerDecomposeLattice(double * vortexlattice, double * vortexremovedlattice, double * fulllattice, double center_omega, double Rtol);

#endif
