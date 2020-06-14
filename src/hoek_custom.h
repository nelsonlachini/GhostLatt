#ifndef HOEK_CUSTOM_H
#define HOEK_CUSTOM_H

double map(double n, int N_ins);

double mapR(double y);

void imposePeriodicity(double * lattice, double * x0, int N_ins);

double getStepBorder( double coord , double dcoord , int N_ins);

double updateStapleBorder(double * lattice , int t, int x, int y, int z, int mi , double * _staple, double * x0, int N_ins);

void coolingInstantonBorderStep(double * lattice, double * x0, int N_ins, int border_length);

double coolInstantonBorder(double * lattice, double * x0, int N_ins, int Ncool, int border_length);

void amuregV(double * amu, double t, double x, double y, double z, int mu, double R, int N_ins);

void creategV(double * g, double t, double x, double y, double z, int N_ins);

void umuregV(double * umu, double t, double x, double y, double z, int mu, int partition_number,double R, int N_ins);

void singtrans(double * uout, double * uin, double t, double x, double y, double z, int mu, int N_ins);

double teperInstanton(double * U, int partition_number, double R, double c0, double * x0, int N_ins, int N_cool);

double teperAntiInstanton(double * U, int partition_number, double R, double c0, double * x0, int N_ins, int N_cool);

#endif
