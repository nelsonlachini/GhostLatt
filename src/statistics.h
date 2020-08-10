#ifndef STATISTICS_H
#define STATISTICS_H

#include "statistics.h"
#include "utilities.h"
#include "algebra.h"

#include "nrutil.h"

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define TOL 1.0e-5

double naive_avg(double * dataSet, int dim);

double naive_sampleVar(double * dataSet, int dim);

int sortInt(int i , int j);

void plaquetteHistogram(double * histo, LatticeLinkSU2 * lattice, int nbins);

void bootstrapSimple(double * output, double * data, int n, int m, int K);

void bootstrapCreutzRatio(double * output, double * wij, double * wiMjM, double * wiMj, double * wijM, int n, int m, int K);

void bootstrapSimpleRatio(double * output, double * wij, double * wijM, int n, int m, int K);

void bootstrapSimpleLog(double * output, double * wij, int n, int m, int K);

void cornellf(double x, double f[], int np);

//from Numerical Recipes 90: functions leading to linear fit()

double gammln(double xx);

void gcf(double *gammcf, double a, double x, double *gln);

void gser(double *gamser, double a, double x, double *gln);

double gammq(double a, double x);

void fit(double x[], double y[], int ndata, double sig[], int mwt, double *a,
    double *b, double *siga, double *sigb, double *chi2, double *q);

void svdvar(double **v, int ma, double w[], double **cvm);

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);

void svdcmp(double **a, int m, int n, double w[], double **v);

double pythag(double a, double b);

void svdfit(double x[], double y[], double sig[], int ndata, double a[],
     int ma, double **u, double **v, double w[], double *chisq,
        void (*funcs)(double, double [], int));
//data blocking
void dataBlockingOnlyDivisors(double * data, int dim);

#endif
