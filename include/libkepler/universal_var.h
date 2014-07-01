#ifndef LIBKEPLER_UNIVERSAL_VAR_H
#define LIBKEPLER_UNIVERSAL_VAR_H

double universal_var_C(double z);
double universal_var_S(double z);
double universal_var_Cseries(double z);
double universal_var_Sseries(double z);
double universal_var_dCdz(double z, double C, double S);
double universal_var_dSdz(double z, double C, double S);
double universal_var_t(
    double sqrtmu,
    double alpha,
    double rv, double r0,
    double x, double z,
    double C, double S);
double universal_var_dtdx(
    double sqrtmu,
    double rv, double r0,
    double x, double z,
    double C, double S);
double universal_var_iterate_x(
    double sqrtmu,
    double alpha,
    double rv, double r0,
    double x0, double t);
double universal_var_guess_x(
    double mu,
    double sqrtmu,
    double alpha,
    double rv, double r0,
    double t);
void universal_var_fg(
    double mu,
    double r0,
    double t,
    double x, double z,
    double C, double S,
    double *f, double *g);
void universal_var_fgdot(
    double sqrtmu,
    double r0,
    double r,
    double x, double z,
    double C, double S,
    double *fdot, double *gdot);

#endif
