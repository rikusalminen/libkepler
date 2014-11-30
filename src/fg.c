#include <libkepler/fg.h>

#include <stdbool.h>
#include <math.h>
#include <float.h>

static bool zero(double x) { return x*x < DBL_EPSILON; }

void fg(
    const double *pos0, const double *vel0,
    const double *pos1, const double *vel1,
    double *f, double *g,
    double *fdot, double *gdot) {
    // specific angular momentum
    double h = pos0[0]*vel0[1] - pos0[1]*vel0[0];

    *f = (pos1[0]*vel0[1] - vel0[0]*pos1[1]) / h;
    *g = (pos0[0]*pos1[1] - pos1[0]*pos0[1]) / h;

    *fdot = (vel1[0]*vel0[1] - vel0[0]*vel1[1]) / h;
    *gdot = (pos0[0]*vel1[1] - vel1[0]*pos0[1]) / h;
}

void fg_series_coefficients(
    double mu,
    double r0,
    double v0,
    double r0v0,
    double *fs,
    double *gs) {
    // f & g series coefficients from
    // Bate, Mueller, White: "Fundamentals of Astrodynamics"
    // See also "Lagrange's fundamental invariants" in
    // Battin: Introduction to the mathematics and methods of astrodynamics
    double u = mu / (r0*r0*r0);
    double p = r0v0 / (r0*r0);
    double q = (v0*v0) / (r0*r0) - u;

    fs[0] = 1.0;
    fs[1] = 0.0;
    fs[2] = -u;
    fs[3] = 3.0*u*p;
    fs[4] = u * (-15.0*p*p + 3.0*q + u);
    fs[5] = 15.0*u*p * (7.0*p*p - 3.0*q - u);
    fs[6] = 0.0; // XXX: more terms makes this less accurate
    /*
    fs[6] = 105.0*u*p*p * (-9.0*p*p + 6.0*q + 2.0*u) -
        u * (45.0*q*q + 24.0*u*p + u*u);
    fs[7] = 315.0*u*p*p*p * (33.0*p*p - 30.0*q - 10*u) +
        63.0*u*p * (25.0*q*q + 14.0*u*p + u*u);
    */

    gs[0] = 0.0;
    gs[1] = 1.0;
    gs[2] = 0.0;
    gs[3] = -u;
    gs[4] = 6.0*u*p;
    gs[5] = u * (-45.0*p*p + 9.0*q + u);
    gs[6] = 30.0*u*p * (14.0*p*p - 6.0*q - u);
    /*
    gs[7] = 315.0*u*p*p * (-15.0*p*p + 10.0*q + 2.0*u) -
        u * (225.0*q*q + 54.0*u*q + u*u);
    gs[8] = 630.0*u*p*p*p * (99.0*p*p - 90.0*q - 20.0*u) +
        126.0*u*p * (75.0*q*q + 24.0*u*q + u*u);
    */
}

void fg_series_t(
    const double *fs,
    const double *gs,
    double dt,
    double *f, double *g,
    double *fdot, double *gdot) {
    // Taylor series using coefficients from fs and gs
    double ff = 0.0, gg = 0.0;
    double ffdot = 0.0, ggdot = 0.0;

    double numer = 1.0, denom = 1.0;
    for(int i = 0; i < FG_SERIES_NUM_COEFFICIENTS; ++i) {
        ffdot += fs[i] * i * numer / denom;
        ggdot += gs[i] * i * numer / denom;

        numer *= i ? dt : 1.0;

        ff += fs[i] * numer / denom;
        gg += gs[i] * numer / denom;

        denom *= i+1;
    }

    *f = ff; *g = gg;
    *fdot = ffdot; *gdot = ggdot;
}

void fg_series(
    double mu,
    double r0,
    double v0,
    double r0v0,
    double dt,
    double *f, double *g,
    double *fdot, double *gdot) {
    double fs[FG_SERIES_NUM_COEFFICIENTS], gs[FG_SERIES_NUM_COEFFICIENTS];

    fg_series_coefficients(mu, r0, v0, r0v0, fs, gs);
    fg_series_t(fs, gs, dt, f, g, fdot, gdot);
}

void fg_eccentric(
    double mu,
    double p, double e,
    double r0,
    double dt, double dE,
    double *f, double *g) {
    double a = p / (1.0 - e*e);

    if(zero(e-1.0)) {
        // parabolic
        double q = p / 2.0;
        *f = 1.0 - (q / r0) * (dE*dE);
        *g = dt - sqrt(p*p*p / mu) * (dE*dE*dE) / 6.0;
    } else if(e < 1.0) {
        // elliptic
        *f = 1.0 - (a/r0) * (1.0 - cos(dE));
        *g = dt - sqrt(a*a*a / mu) * (dE - sin(dE));
    } else if(e > 1.0) {
        // hyperbolic
        *f = 1.0 - (a/r0) * (1.0 - cosh(dE));
        *g = dt - sqrt(-a*a*a / mu) * (sinh(dE) - dE);
    }
}

void fgdot_eccentric(
    double mu,
    double p, double e,
    double r0, double r,
    double dE,
    double *fdot, double *gdot) {
    double a = p / (1.0 - e*e);

    if(zero(e-1.0)) {
        // parabolic
        double q = p / 2.0;
        *fdot = -sqrt(mu * p) * dE / (r0 * r);
        *gdot = 1.0 - (q / r) * dE*dE;
    } else if(e < 1.0) {
        // elliptic
        *fdot = - sqrt(mu * a) / (r0 * r) * sin(dE);
        *gdot = 1.0 - (a / r) * (1.0 - cos(dE));
    } else if(e > 1.0) {
        // hyperbolic
        *fdot = - sqrt(mu * -a) / (r0 * r) * sinh(dE);
        *gdot = 1.0 - (a / r) * (1.0 - cosh(dE));
    }
}

void fg_true(
    double mu,
    double p,
    double r0, double r,
    double df,
    double *f, double *g) {
    *f = 1.0 - r/p * (1.0 - cos(df));
    *g = r*r0 * sin(df) / sqrt(mu * p);
}

void fgdot_true(
    double mu,
    double p,
    double r0, double r,
    double df,
    double *fdot, double *gdot) {
    *fdot = sqrt(mu / p) * tan(df/2.0) * ((1.0 - cos(df))/p - 1.0/r - 1.0/r0);
    *gdot = 1.0 - r0/p * (1.0 - cos(df));
}
