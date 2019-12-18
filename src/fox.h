#ifndef FOX_H
#define FOX_H

double xymapFox(double y , double x0);

double xymapRadius(double y , double x0);

int isoutside(double t,double x,double y,double z, double x0, int M);

void amuregVFox(double * amu, double t, double x, double y, double z, int mu, double x0,double R);

void amusingVFox(double * amu, double t , double x, double y, double z, int mi, double x0, double R);

void umuregVFox(double * umu, int t, int x, int y, int z, int mi, double x0, double R);

void umusingVFox(double * umu, int t, int x, int y, int z, int mi, double x0, double R);

void creategFoxV(double * g, double t, double x, double y, double z, double x0);

void foxInstanton(double * U , double R, double x0 , int M);

#endif
