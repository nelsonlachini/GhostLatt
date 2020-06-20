#ifndef INVERTERS_H
#define INVERTERS_H

#include <complex.h>

#include "global.h"

// extern const int ITERATION_TOL;

//CG and PM with real vectors
double rfpapply(double * Mvout , LatticeSU2 * U , double * vin , int dim);

double rfpapplyOrthogonal(double * Mvout , LatticeSU2 * lattice , double * vin , int dim);

double rfpcgmr(double * x , LatticeSU2 * lattice , double * b,  int dim, double tol, int * iCG, double * r_CG);

void rgetOrthogonalSpace(double * c, double * b, LatticeSU2 * lattice, double cg_tol, int * iCG, double * r_CG);

int rsmallest_eigen_cg(LatticeSU2 * lattice, double * lambda1, double * eigen_out, double eigen_tol 
    , double * vguess, double * r_cg, int * iCG, double cg_tol);

//CG and PM with complex vector
double fpapply(double complex * Mvout , LatticeSU2 * lattice , double complex * vin , int dim);

double fpcgmr(double complex * x , LatticeSU2 * lattice , double complex * b,  int dim);

void getOrthogonalSpace(double complex * c, double complex * b, LatticeSU2 * lattice);

double complex smallest_eigen_cg(LatticeSU2 * lattice, double * lambda1, double complex * eigen_out);

//biCGstab with complex vectors
void fpbicgstabmr(double complex * x , LatticeSU2 * lattice , double complex * b,  int dim, double cg_tol, double * r_cg, int * icg);

void bigetOrthogonalSpace(double complex * c, double complex * b, LatticeSU2 * lattice, double cg_tol, double * r_cg, int * icg);

double smallest_eigen_bicgstab(LatticeSU2 * lattice, double * lambda1, double * eigen_out, double eigen_tol, double * r_cg, int * icg, double cg_tol);

//R functional computation
double  FsecondDerivative(LatticeSU2 * lattice, double complex * eigen_vector_out);

double rFThirdFourthDerivative(LatticeSU2 * lattice, double * ev_in, double * third, double * fourth);

//OTHER STUFF
double complex largest_eigen_cg(LatticeSU2 * lattice, double * lambda1, double complex * eigen_out);

void divergentConfigProcedure(LatticeSU2 * escaled_lattice, double new_delta_scale, double * tau, int * divergenConfigN, double r_CG, int iCG, double * lambda1, double * ev_guess, double eigen_tol, double cg_tol);

double horizonWalk(LatticeSU2 * target_lattice, double * nHorizonOut, double * rhoHorizonOut, double * ev_guess
		, double * lambda1Out, double * rAfterOut, double * rBeforeOut, double * thirdOut, double * thirdAbsOut, double * fourthOut
		, double * rOut, double * pwProjOut, double eigen_tol, int * divergentConfig, double cg_tol);

double verifyOrthogonalization(LatticeSU2 * lattice, double * vin);
 #endif
