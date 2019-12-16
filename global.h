#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "algebra.h"
#include "utilities.h"

#define MAXIMUM_NAME_LENGTH 80
#define NBIN 50

extern const int N;
extern const int Nt;
extern const int N2;
extern const int spatialV;
extern const int totalV;
extern const int colorV;
extern const int dimLattice;

extern const double beta;

extern long * global_seed;

#endif