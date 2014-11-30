#include <libkepler/kepler.h>
#include <libkepler/fg.h>

#include <math.h>
#include <float.h>

#include "numtest.h"

static double square(double x) { return x*x; }
static double cube(double x) { return x*x*x; }

static double dot(const double *a, const double *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double mag(const double *a) { return sqrt(dot(a, a)); }

static bool eqv3(double *a, double *b) {
    double diff[3] = { a[0]-b[0], a[1]-b[1], a[2]-b[2] };
    return (ZEROF(dot(a,a)) && ZEROF(dot(b,b))) ||
        ZEROF(dot(diff, diff)/(dot(a, a) + dot(b, b)));
}

void fg_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 5, "");

    double p = 1.0 + params[0] * 1.0e8;
    double e = params[1] * 4.0;
    double mu = 1.0 + params[2] * 1.0e8;
    double a = p / (1.0 - square(e));
    double n = ZEROF(e - 1.0) ?
        sqrt(mu / cube(p)) :    // parabolic
        sqrt(mu / cube(fabs(a)));
    double t0 = 0.0;

    struct kepler_elements elements = { p, e, n, 0.0, 0.0, 0.0, t0 };

    double maxf = e > 1.0 ?
        acos(1.0/e) :       // hyperbolic
        (ZEROF(e-1.0) ?
            4*M_PI/5 :      // parabolic
            M_PI);          // closed orbit

    double k1 = (-1.0 + params[3] * 2.0);
    double f1 = maxf * k1;
    double E1 = kepler_anomaly_true_to_eccentric(e, f1);
    double M1 = kepler_anomaly_eccentric_to_mean(e, E1);
    double t1 = t0 + M1 / n;
    double r1 = p / (1.0 + e * cos(f1));

    double k2 = (-1.0 + params[4] * 2.0);
    double f2 = maxf * k2;
    double E2 = kepler_anomaly_true_to_eccentric(e, f2);
    double M2 = kepler_anomaly_eccentric_to_mean(e, E2);
    double t2 = t0 + M2 / n;
    double r2 = p / (1.0 + e * cos(f2));

    double pos1[3], vel1[3], pos2[3], vel2[3];
    kepler_elements_to_state_E(&elements, E1, pos1, vel1);
    kepler_elements_to_state_E(&elements, E2, pos2, vel2);

    double f, g, fdot, gdot;
    fg(pos1, vel1, pos2, vel2, &f, &g, &fdot, &gdot);

    ASSERT(isfinite(f) && isfinite(g) && isfinite(fdot) && isfinite(gdot),
        "fg not NaN");

    double pos[3], vel[3];
    for(int i = 0; i < 3; ++i) {
        pos[i] = f*pos1[i] + g*vel1[i];
        vel[i] = fdot*pos1[i] + gdot*vel1[i];
    }

    ASSERT(eqv3(pos, pos2), "fg position identity");
    ASSERT(eqv3(vel, vel2), "fg velocity identity");

    ASSERT_EQF(f*gdot, 1.0 + fdot*g, "fg identity (sanity)");


    // Eccentric anomaly
    double fE, gE;
    fg_eccentric(mu, p, e, r1, t2-t1, E2-E1, &fE, &gE);

    ASSERT(isfinite(fE) && isfinite(gE),
        "fg not NaN (eccentric anomaly)");

    double posE[3];
    for(int i = 0; i < 3; ++i)
        posE[i] = fE*pos1[i] + gE*vel1[i];
    ASSERT(eqv3(posE, pos2), "fg position (eccentric anomaly)");

    double fdotE, gdotE;
    double rE = mag(posE);
    fgdot_eccentric(mu, p, e, r1, rE, E2-E1, &fdotE, &gdotE);

    ASSERT(isfinite(fdotE) && isfinite(gdotE),
        "fgdot not NaN (eccentric anomaly)");

    double velE[3];
    for(int i = 0; i < 3; ++i)
        velE[i] = fdotE*pos1[i] + gdotE*vel1[i];
    ASSERT(eqv3(velE, vel2), "fg velocity (eccentric anomaly)");

    ASSERT_EQF(fE*gdotE, 1.0 + fdotE*gE, "fg identity (eccentric anomaly)");

    // True anomaly
    double ff, gf;
    fg_true(mu, p, r1, r2, f2-f1, &ff, &gf);

    ASSERT(isfinite(ff) && isfinite(gf),
        "fg not NaN (true anomaly)");

    double posf[3];
    for(int i = 0; i < 3; ++i)
        posf[i] = ff*pos1[i] + gf*vel1[i];
    ASSERT(eqv3(posf, pos2), "fg position (true anomaly)");

    double fdotf, gdotf;
    fgdot_true(mu, p, r1, r2, f2-f1, &fdotf, &gdotf);

    ASSERT(isfinite(fdotf) && isfinite(gdotf),
        "fgdot not NaN (true anomaly)");

    double velf[3];
    for(int i = 0; i < 3; ++i)
        velf[i] = fdotf*pos1[i] + gdotf*vel1[i];
    if(fabs(f2-f1) < M_PI*7.0/8.0) // inaccurate near 180 degrees
        ASSERT(eqv3(velf, vel2), "fg velocity (true anomaly)");

    ASSERT_EQF(ff*gdotf, 1.0 + fdotf*gf, "fg identity (eccentric anomaly)");

    // Time series
    double v1 = mag(vel1), r1v1 = dot(pos1, vel1);
    double fs, gs, fdots, gdots;
    double f_coeffs[FG_SERIES_NUM_COEFFICIENTS], g_coeffs[FG_SERIES_NUM_COEFFICIENTS];

    fg_series_coefficients(mu, r1, v1, r1v1, f_coeffs, g_coeffs);
    for(int i = 0; i < FG_SERIES_NUM_COEFFICIENTS; ++i)
        ASSERT(isfinite(f_coeffs[i]) && isfinite(g_coeffs[i]),
            "fg time series coefficients not NaN");
    fg_series_t(f_coeffs, g_coeffs, t2-t1, &fs, &gs, &fdots, &gdots);

    ASSERT(isfinite(fs) && isfinite(gs) && isfinite(fdots) && isfinite(gdots),
        "fg not NaN (time series)");

    double poss[3], vels[3];
    for(int i = 0; i < 3; ++i) {
        poss[i] = fs*pos1[i] + gs*vel1[i];
        vels[i] = fdots*pos1[i] + gdots*vel1[i];
    }

    // time series is only accurate for small angles
    ASSERT(fabs(f2-f1) > 5*M_PI/180.0 ||
        eqv3(poss, pos2), "fg position (time series)");
    ASSERT(fabs(f2-f1) > 3*M_PI/180.0 ||
        eqv3(vels, vel2), "fg velocity (time series)");

    ASSERT(fabs(f2-f1) > 5*M_PI/180.0 ||
        EQF(fs*gdots, 1.0 + fdots*gs), "fg identity (time series)");
}
