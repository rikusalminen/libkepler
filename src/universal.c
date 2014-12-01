#include <libkepler/universal.h>

double universal_r_cs(
    double mu,
    double alpha,
    double s,
    double r0, double r0v0,
    const double *cs) {
    (void)alpha;
    return r0 * cs[0] + r0v0 * s * cs[1] + mu * s*s * cs[2];
}

double universal_t_cs(
    double mu,
    double alpha,
    double s,
    double r0, double r0v0,
    const double *cs) {
    (void)alpha;
    return r0 * s * cs[1] + r0v0 * s*s * cs[2] + mu * s*s*s * cs[3];
}

void universal_fg_cs(
    double mu,
    double alpha,
    double s,
    double r0, double r0v0,
    const double *cs,
    double *f, double *g) {
    (void)alpha;
    *f = 1.0 - (mu/r0) * s*s * cs[2];
    *g = r0 * s * cs[1] + r0v0 * s*s * cs[2];
    // *g = dt - mu * s*s*s * cs[3];
}


void universal_fgdot_cs(
    double mu,
    double alpha,
    double s,
    double r0, double r,
    const double *cs,
    double *fdot, double *gdot) {
    (void)alpha;
    *fdot = -mu / (r0*r) * s * cs[1];
    *gdot = 1.0 - (mu/r) * s*s * cs[2];
}
