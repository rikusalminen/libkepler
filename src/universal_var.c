#include <stdbool.h>
#include <float.h>
#include <math.h>

static double dot(const double *a, const double *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double mag(const double *a) {
    return sqrt(dot(a, a));
}

static bool zero(double x) {
    return x*x < DBL_EPSILON;
}

static double sign(double x) { return x < 0 ? -1.0 : 1.0; }
static double square(double x) { return x*x; }
static double cube(double x) { return x*x*x; }


double universal_var_C(double z) {
    if(zero(z))
        return 1.0/2.0;
    if(z < 0.0)
        return (1.0 - cosh(sqrt(-z)))/z;
    return (1.0 - cos(sqrt(z)))/z;
}

double universal_var_S(double z) {
    if(zero(z))
        return 1.0/6.0;
    if(z < 0.0) {
        double root = sqrt(-z);
        return (sinh(root) - root)/sqrt(-z*z*z);
    }

    double root = sqrt(z);
    return (root - sin(root)) / sqrt(z*z*z);
}

double universal_var_Cseries(double z) {
    const int num_steps = 30;
    double threshold = DBL_EPSILON;

    int i = 1;
    double k = 1.0/2.0, c = k;
    do {
        k *= -z / ((2*i+1) * (2*i+2));
        c += k;
    } while(++i <= num_steps && k*k > threshold);

    return c;
}

double universal_var_Sseries(double z) {
    const int num_steps = 30;
    double threshold = DBL_EPSILON;

    int i = 1;
    double k = 1.0/6.0, s = k;
    do {
        k *= -z / ((2*i+2)*(2*i+3));
        s += k;
    } while(++i <= num_steps && k*k > threshold);

    return s;
}

double universal_var_dCdz(double z, double C, double S) {
    return (1.0 - z*S - 2.0*C)/(2.0*z);
}

double universal_var_dSdz(double z, double C, double S) {
    return (C - 3.0*S)/(2.0*z);
}

double universal_var_t(
    double sqrtmu,
    double alpha,
    double rv, double r0,
    double x, double z,
    double C, double S) {
    (void)z;
    return (rv/sqrtmu * square(x) * C + (1.0 - r0*alpha) * cube(x) * S + r0*x)/sqrtmu;
}

double universal_var_dtdx(
    double sqrtmu,
    double rv, double r0,
    double x, double z,
    double C, double S) {
    return (square(x) * C + (rv/sqrtmu)*(1.0 - z*S) + r0 * (1.0 - z*C))/sqrtmu;
}

double universal_var_iterate_x(
    double sqrtmu,
    double alpha,
    double rv, double r0,
    double x0, double t) {
    const int num_steps = 30;
    double threshold = DBL_EPSILON;

    int i = 1;
    double k = 0.0/0.0;
    double x = x0;
    do {
        double z = alpha * square(x);
        double C = universal_var_Cseries(z);
        double S = universal_var_Sseries(z);
        double tt = universal_var_t(sqrtmu, alpha, rv, r0, x, z, C, S);
        double dtdx = universal_var_dtdx(sqrtmu, rv, r0, x, z, C, S);

        k = (t - tt) / dtdx;
        x += k;
    } while(++i <= num_steps && k*k > threshold);

    return x;
}

double universal_var_guess_x(
    double mu,
    double sqrtmu,
    double alpha,
    double rv, double r0,
    double t) {
    //if(zero(alpha))
        //return 0.0 / 0.0; // XXX: solve X from parabolic anomaly?!

    if(alpha < 0.0) {   // hyperbolic trajectory
        double sqrta = sqrt(-1.0 / alpha);

        return sign(t)*sqrta *
            log((-2.0*mu*t*alpha)/(rv + sign(t)*sqrtmu*sqrta*(1.0 - r0*alpha)));
    }

    return sqrtmu * alpha * t;
}

void universal_var_fg(
    double sqrtmu,
    double r0,
    double t,
    double x, double z,
    double C, double S,
    double *f, double *g) {
    (void)z;

    *f = 1.0 - square(x)/r0 * C;
    *g = t - cube(x)/sqrtmu * S;
}

void universal_var_fgdot(
    double sqrtmu,
    double r0,
    double r,
    double x, double z,
    double C, double S,
    double *fdot, double *gdot) {
    *fdot = sqrtmu/(r0*r) * x * (z*S - 1.0);
    *gdot = 1.0 - square(x)/r * C;
}

void universal_var_fg_E(
    double sqrtmu,
    double a,
    double e,
    double r0,
    double delta_t, double delta_E,
    double *f, double *g) {
    // XXX: parabolic trajectory!
    if(e > 1.0) {
        *f = 1.0 - a/r0 * (1.0 - cos(delta_E));
        *g = delta_t - (pow(-a, 3.0/2.0) / sqrtmu) * (sinh(delta_E) - delta_E);
    } else {
        *f = 1.0 - a/r0 * (1.0 - cos(delta_E));
        *g = delta_t - (pow(a, 3.0/2.0) / sqrtmu) * (delta_E - sin(delta_E));
    }
}

void universal_var_fgdot_E(
    double sqrtmu,
    double a,
    double e,
    double r0,
    double r,
    double delta_t, double delta_E,
    double *fdot, double *gdot) {
    (void)delta_t;
    // XXX: parabolic trajectory!
    if(e > 1.0) {
        *fdot = -sqrtmu*sqrt(-a)*sinh(delta_E) / (r*r0);
        *gdot = 1.0 - a/r * (1.0 - cosh(delta_E));
    } else {
        *fdot = -sqrtmu*sqrt(a)*sin(delta_E)/(r*r0);
        *gdot = 1.0 - a/r * (1.0 - cos(delta_E));
    }
}

void universal_var_fg_f(
    double sqrtmu,
    double p,
    double r0, double r,
    double delta_f,
    double *f, double *g) {
    *f = 1.0 - r/p * (1.0 - cos(delta_f));
    *g = r*r0 * sin(delta_f) / (sqrtmu * sqrt(p));
}

void universal_var_fgdot_f(
    double sqrtmu,
    double p,
    double r0, double r,
    double delta_f,
    double *fdot, double *gdot) {
    *fdot = sqrtmu/sqrt(p) * tan(delta_f/2.0) * ((1.0 - cos(delta_f))/p - 1.0/r - 1.0/r0);
    *gdot = 1.0 - r0/p * (1.0 - cos(delta_f));
}

void universal_var(
    double mu,
    const double *pos0, const double *vel0,
    double t,
    double *pos, double *vel) {
    double r0 = mag(pos0), v02 = dot(vel0, vel0); //, v0 = sqrt(v02);
    double rv = dot(pos0, vel0);
    double sqrtmu = sqrt(mu);
    double alpha = (2.0*mu/r0 - v02) / mu;

    double x0 = universal_var_guess_x(mu, sqrtmu, alpha, rv, r0, t);
    double x = universal_var_iterate_x(sqrtmu, alpha, rv, r0, x0, t);
    double z = alpha * square(x);

    double C = universal_var_Cseries(z);
    double S = universal_var_Sseries(z);

    double f, g;
    universal_var_fg(sqrtmu, r0, t, x, z, C, S, &f, &g);

    for(int i = 0; i < 3; ++i)
        pos[i] = f*pos0[i] + g*vel0[i];
    double r = mag(pos);

    double fdot, gdot;
    universal_var_fgdot(sqrtmu, r0, r, x, z, C, S, &fdot, &gdot);

    for(int i = 0; i < 3; ++i)
        vel[i] = fdot*pos0[i] + gdot*vel0[i];
}
