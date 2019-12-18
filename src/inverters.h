#ifndef INVERTERS_H
#define INVERTERS_H

#include <complex.h>

// extern const int ITERATION_TOL;

//CG and PM with real vectors
double rfpapply(double * Mvout , double * U , double * vin , int dim);

double rfpapplyOrthogonal(double * Mvout , double * U , double * vin , int dim);

double rfpcgmr(double * x , double * lattice , double * b,  int dim, double tol, int * iCG, double * r_CG);

void rgetOrthogonalSpace(double * c, double * b, double * lattice, double cg_tol, int * iCG, double * r_CG);

int rsmallest_eigen_cg(double * lattice, double * lambda1, double * eigen_out, double eigen_tol 
    , double * vguess, double * r_cg, int * iCG, double cg_tol);

//CG and PM with complex vector
double fpapply(double complex * Mvout , double * U , double complex * vin , int dim);

double fpcgmr(double complex * x , double * lattice , double complex * b,  int dim);

void getOrthogonalSpace(double complex * c, double complex * b, double * lattice);

double complex smallest_eigen_cg(double * lattice, double * lambda1, double complex * eigen_out);

//biCGstab with complex vectors
void fpbicgstabmr(double complex * x , double * lattice , double complex * b,  int dim, double cg_tol, double * r_cg, int * icg);

void bigetOrthogonalSpace(double complex * c, double complex * b, double * lattice, double cg_tol, double * r_cg, int * icg);

double smallest_eigen_bicgstab(double * lattice, double * lambda1, double * eigen_out, double eigen_tol, double * r_cg, int * icg, double cg_tol);

//R functional computation
double  FsecondDerivative(double * lattice, double complex * eigen_vector_out);

double rFThirdFourthDerivative(double * lattice, double * ev_in, double * third, double * fourth);

//OTHER STUFF
double complex largest_eigen_cg(double * lattice, double * lambda1, double complex * eigen_out);

void divergentConfigProcedure(double * escaled_lattice, double new_delta_scale, double * tau, int * divergenConfigN, double r_CG, int iCG, double * lambda1, double * ev_guess, double eigen_tol, double cg_tol);

double horizonWalk(double * target_lattice, double * nHorizonOut, double * rhoHorizonOut, double * ev_guess
		, double * lambda1Out, double * rAfterOut, double * rBeforeOut, double * thirdOut, double * thirdAbsOut, double * fourthOut
		, double * rOut, double * pwProjOut, double eigen_tol, int * divergentConfig, double cg_tol);

double verifyOrthogonalization(double * lattice, double * vin);
 #endif
