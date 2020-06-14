#ifndef GAUGE_FIX_H
#define GAUGE_FIX_H


//####################################################IN TEST
void reunitg(double * g);

double calcEps(double * U , double * g);

double calcFunctional(double * U);

double * getLatticeReal(double * real , int t , int x , int y , int z);

double latticeAvgReal(double * var);

//////////////////////////auxiliars
void initg(double * g);

double applyFix(double * Uout, double * Uin , double * g);

void copyg(double * gout, double * gin);

double su2distance(double * U, double * W);

//SERIAL VERSIONS

double los_alamos(double * U , double * g , double * e2);

double cornell(double * U , double * g, double * e2 , double alpha);

double overrelaxation(double * U , double * g, double * e2 , double omega);

double stochastic(double * U , double * g , double * e2, double p);

// double fourier(double * U , double * g , double * e2, double alpha);

//AUTOMATIZATION
double fixLatticeStoch(double * lout, double * lin, double * g, double p_stoch, double e2tol);

double fixLatticeLosalamos(double * lout, double * lin, double * g , double e2tol);

double fixLatticeOver(double * lout, double * lin, double * g, double omega, double e2tol);

//AUTO CALIBRATION

double calce2(double * U);

double calibrate_over(double * lattice , double tol);

double calibrate_cornell(double * lattice , double tol, double step_size, double lower_alpha, double upper_alpha);

double calibrate_stoc(double * lattice , double tol);

double * getQ(int ni , int xni);

double * getAvgQ(int ni);

void calcAvgQ(int ni , double * sum);

void calcQni( int ni  , double * lattice);

//double calce6(double * lattice , double * g);

#endif
