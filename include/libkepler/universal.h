#ifndef LIBKEPLER_UNIVERSAL_H
#define LIBKEPLER_UNIVERSAL_H

double universal_r_cs(
    double mu,
    double alpha,
    double s,
    double r0, double r0v0,
    const double *cs);

double universal_t_cs(
    double mu,
    double alpha,
    double s,
    double r0, double r0v0,
    const double *cs);

void universal_fg_cs(
    double mu,
    double alpha,
    double s,
    double r0, double r0v0,
    const double *cs,
    double *f, double *g);

void universal_fgdot_cs(
    double mu,
    double alpha,
    double s,
    double r0, double r0v0,
    const double *cs,
    double *fdot, double *gdot);

double universal_guess_s();
double universal_iterate_s();

#endif
