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
    double sqrtmu,
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
void universal_var_fg_E(
    double sqrtmu,
    double a,
    double e,
    double r0,
    double delta_t, double delta_E,
    double *f, double *g);
void universal_var_fgdot_E(
    double sqrtmu,
    double a,
    double e,
    double r0,
    double r,
    double delta_t, double delta_E,
    double *fdot, double *gdot);
void universal_var_fg_f(
    double sqrtmu,
    double p,
    double r0, double r,
    double delta_f,
    double *f, double *g);
void universal_var_fgdot_f(
    double sqrtmu,
    double p,
    double r0, double r,
    double delta_f,
    double *fdot, double *gdot);

void universal_var(
    double mu,
    const double *pos0, const double *vel0,
    double t,
    double *pos, double *vel);

double universal_var_stumpff_c0(double alpha, double s);
double universal_var_stumpff_c1(double alpha, double s);
double universal_var_stumpff_c2(double alpha, double s);
double universal_var_stumpff_c3(double alpha, double s);
double universal_var_stumpff_series(int k, double z);

void universal_var_stumpff_fast(double z, double *cs);

#endif
