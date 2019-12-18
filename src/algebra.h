#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <complex.h>

int getSignal(double number);

////////////////////////////// SU2 MATRIX/VECTOR ALGEBRA

double tracev( double * pvector );

void sumv( double * pout , double * B , double * C );

void mmprodv( double * pout , double  * B, double  * C );


void cmprodv( double * pout , double k, double * C );   //needs revising

void hermcv(double * Adagger , double * A);

double detv(double * A);

void setidv(double * vector);

void matrixexpv(double * vout, double * vin, int order);

void setzerov(double * vector);

void su2sroot(double * out, double * in);

////////////////////////////// GENERAL REAL MATRIX ALGEBRA

double tracer( double * pmatrix , int dim);

void setidmr(double * M , int dim);

void setzeromr(double * M , int dim);

void sumr( double * pout , double * B , double * C , int dim);

void mmprodr( double * mout , double * B , double * C , int dim);

void transpmr(double * Atransposed , double * A , int dim);

void setzerovr(double * vr , int dim);

double inprodvr(double * a , double * b , int dim);

double normvr(double * a , int dim);

void cprodvr(double * vout , double k , double * v , int dim);

void sumvr(double * vout , double * v1 , double * v2 , int dim);

void reunitvr(double * vector, int dim);

void mvrprod(double * vout , double * B , double * vcolumn , int dim);


////////////////////////////// GENERAL COMPLEX MATRIX ALGEBRA

void setzerovc(double complex * vc , int dim);

double complex tracec( double complex * pmatrix , int dim);

void setidmc(double complex * M , int dim);

void setzeromc(double complex * M , int dim);

void sumc( double complex * mout , double complex * B , double complex * C , int dim);

void mmprodc( double complex * Mout , double complex * B , double complex * C , int dim);

void cmprodc(double complex * mout , double complex k , double complex * C , int dim);

void mvprodr(double complex * vout , double * B , double complex * vcolumn , int dim);

void vmprodr(double complex * vout , double complex* vline , double * B , int dim);

void sumvc(double complex * vout , double complex * v1 , double complex * v2 , int dim);

void cprodvc(double complex * vout , double complex k , double complex * v , int dim);

void hermcvc(double complex * Vdagger , double complex * V , int dim);

void hermcmc(double complex * Mdagger , double complex * M , int dim);

double complex inprodvc(double complex * a , double complex * b , int dim);

double normvc(double complex * a , int dim);

void reunitvc(double complex * vector, int dim);

void mvprodc(double complex * vout , double complex * B , double complex * vcolumn , int dim);

void vmprodc(double complex * vout , double complex * vline , double complex * B , int dim);

double complex det2mc(double complex * v);


/////////////////////////////// LINEAR SYSTEMS ALGORITHMS

void bicgmr(double complex * x , double * A , double complex * b , int dim);

void bicgmc(double complex * x , double complex * A , double complex * b , int dim);

void cgmr(double * x , double * A , double * b , int dim);

void cgmc(double complex * x , double complex * A , double complex * b , int dim);


///////////////////////////////

void setzeroqvector(int * v);

int aseps(int i , int j , int k);

//////////////////////////////////////NOT INCORPORATED YET

void setonevc(double complex * vc,int dim);

void setrandomvc(double complex * vc,int dim);

void setonevr(double * vr,int dim);

void setplanewave(double complex * target, double * p , int atarget);

void setplanewaveallcolors(double complex * target, double * p );

void rsetplanewave(double * target, double * p , int atarget);

void rsetplanewaveallcolors(double * target, double * p);

void setrandomvr(double * vr,int dim);

void invMv(double * vout, double * vin);

void invMc(double complex * vout, double complex * vin);

int eps4(int a, int b , int c , int d);

int eps4tilde(int a, int b, int c, int d);

#endif
