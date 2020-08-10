#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>

#define MAXIMUM_NAME_LENGTH 80
#define NBIN 50

extern const int N;

extern const int Nt;
extern const int N2;
extern const int spatialV;
extern const int totalV;
extern const int colorV;
extern const int dimLattice;

struct LatticeLinkSU2{
    int N;
    int Nt;
    int N2;
    int spatialV;
    int totalV;
    int colorV;
    int dimLattice;
    double * U;
};

struct LatticeColorVectorReal{
    int N;
    int Nt;
    int dimVector;
    double * vec;
};

struct LatticeColorVectorComplex{
    int N;
    int Nt;
    int dimVector;
    double complex * vec;
};

struct LatticeGaugeSU2{
    int N;
    int N2;
    int Nt;
    int spatialV;
    int totalV;
    int dimg;
    double * g;
};

typedef struct LatticeLinkSU2 LatticeLinkSU2;
typedef struct LatticeColorVectorReal LatticeColorVectorReal;
typedef struct LatticeColorVectorComplex LatticeColorVectorComplex;
typedef struct LatticeGaugeSU2 LatticeGaugeSU2;

extern long * global_seed;

#endif
