#ifndef HOEK_CUSTOM_H
#define HOEK_CUSTOM_H

double map(double n, int N_ins);

double mapR(double y);

void imposePeriodicity(double * lattice, double * x0, int N_ins);

double coolInstantonBorder(double * lattice, double * x0, int N_ins, int Ncool, int border_length);

double teperInstanton(double * U, int partition_number, double R, double c0, double * x0, int N_ins, int N_cool);

double teperAntiInstanton(double * U, int partition_number, double R, double c0, double * x0, int N_ins, int N_cool);

#endif
