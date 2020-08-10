#ifndef INVERTERS_H
#define INVERTERS_H

#include <complex.h>

#include "global.h"

// extern const int ITERATION_TOL;

//CG and PM with real vectors
double rfpapply(LatticeColorVectorReal * Mvout , LatticeLinkSU2 * lattice , LatticeColorVectorReal * vin);

double rfpapplyOrthogonal(LatticeColorVectorReal * Mvout , LatticeLinkSU2 * lattice , LatticeColorVectorReal * vin);

double rfpcgmr(LatticeColorVectorReal * x , LatticeLinkSU2 * lattice , LatticeColorVectorReal * b,  int dim, double tol, int * iCG, double * r_CG);

void rgetOrthogonalSpace(LatticeColorVectorReal * c, LatticeColorVectorReal * b, LatticeLinkSU2 * lattice, double cg_tol, int * iCG, double * r_CG);

int rsmallest_eigen_cg(LatticeLinkSU2 * lattice, double * lambda1, LatticeColorVectorReal * eigen_out, double eigen_tol 
    , LatticeColorVectorReal * vguess, double * r_cg, int * iCG, double cg_tol);

//CG and PM with complex vector
double fpapply(LatticeColorVectorComplex * Mvout , LatticeLinkSU2 * lattice , LatticeColorVectorComplex * vin);

double fpcgmr(LatticeColorVectorComplex * x , LatticeLinkSU2 * lattice , LatticeColorVectorComplex * b);

void getOrthogonalSpace(LatticeColorVectorComplex * c, LatticeColorVectorComplex * b, LatticeLinkSU2 * lattice);

double complex smallest_eigen_cg(LatticeLinkSU2 * lattice, double * lambda1, LatticeColorVectorComplex * eigen_out);

//biCGstab with complex vectors
void fpbicgstabmr(LatticeColorVectorComplex * x , LatticeLinkSU2 * lattice , LatticeColorVectorComplex * b, double cg_tol, double * r_cg, int * icg);

void bigetOrthogonalSpace(LatticeColorVectorComplex * c, LatticeColorVectorComplex * b, LatticeLinkSU2 * lattice, double cg_tol, double * r_cg, int * icg);

// double smallest_eigen_bicgstab(LatticeLinkSU2 * lattice, double * lambda1, double * eigen_out, double eigen_tol, double * r_cg, int * icg, double cg_tol);

//R functional computation
double  FsecondDerivative(LatticeLinkSU2 * lattice, LatticeColorVectorComplex * eigen_vector_out);

double rFThirdFourthDerivative(LatticeLinkSU2 * lattice, double * ev_in, double * third, double * fourth);

//OTHER STUFF
double complex largest_eigen_cg(LatticeLinkSU2 * lattice, double * lambda1, LatticeColorVectorComplex * eigen_out);

void divergentConfigProcedure(LatticeLinkSU2 * escaled_lattice, double new_delta_scale, double * tau, int * divergenConfigN, double r_CG, int iCG, double * lambda1, LatticeColorVectorReal * ev_guess, double eigen_tol, double cg_tol);

double horizonWalk(LatticeLinkSU2 * target_lattice, double * nHorizonOut, double * rhoHorizonOut, double * ev_guess
		, double * lambda1Out, double * rAfterOut, double * rBeforeOut, double * thirdOut, double * thirdAbsOut, double * fourthOut
		, double * rOut, double * pwProjOut, double eigen_tol, int * divergentConfig, double cg_tol);

double verifyOrthogonalization(LatticeLinkSU2 * lattice, LatticeColorVectorReal * vin);

 #endif
