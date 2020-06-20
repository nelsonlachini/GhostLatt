#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define MAXIMUM_NAME_LENGTH 80
#define NBIN 50

extern const int N;

extern const int Nt;
extern const int N2;
extern const int spatialV;
extern const int totalV;
extern const int colorV;
extern const int dimLattice;

struct LatticeSU2{
    int N;
    int Nt;
    int N2;
    int spatialV;
    int totalV;
    int colorV;
    int dimLattice;
    double * U;
};

typedef struct LatticeSU2 LatticeSU2;

extern long * global_seed;

#endif
