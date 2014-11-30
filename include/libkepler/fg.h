#ifndef LIBKEPLER_FG_H
#define LIBKEPLER_FG_H

void fg(
    const double *pos0, const double *vel0,
    const double *pos1, const double *vel1,
    double *f, double *g,
    double *fdot, double *gdot);

#define FG_SERIES_NUM_COEFFICIENTS 7

void fg_series_coefficients(
    double mu,
    double r0,
    double v0,
    double r0v0,
    double *fs,
    double *gs);
void fg_series_t(
    const double *fs,
    const double *gs,
    double dt,
    double *f, double *g,
    double *fdot, double *gdot);
void fg_series(
    double mu,
    double r0,
    double v0,
    double r0v0,
    double dt,
    double *f, double *g,
    double *fdot, double *gdot);

void fg_eccentric(
    double mu,
    double p, double e,
    double r0,
    double dt, double dE,
    double *f, double *g);
void fgdot_eccentric(
    double mu,
    double p, double e,
    double r0, double r,
    double dE,
    double *fdot, double *gdot);

void fg_true(
    double mu,
    double p,
    double r0, double r,
    double df,
    double *f, double *g);
void fgdot_true(
    double mu,
    double p,
    double r0, double r,
    double df,
    double *f, double *g);

#endif
