#include <math.h>

#include <libkepler/kepler.h>
#include <libkepler/universal_var.h>

#include "numtest.h"

static double dot(const double *a, const double *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double mag(const double *a) {
    return sqrt(dot(a, a));
}

static double square(double x) { return x*x; }
static double cube(double x) { return x*x*x; }

static bool eqv3(double *a, double *b) {
    double diff[3] = { a[0]-b[0], a[1]-b[1], a[2]-b[2] };
    return (ZEROF(dot(a,a)) && ZEROF(dot(b,b))) ||
        ZEROF(dot(diff, diff)/(dot(a, a) + dot(b, b)));
}

void universal_var_stumpff_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 2, "");

    const double maxs = 4.0*M_PI, maxalpha = 1.0;
    double s = (-1.0 + 2.0 * params[0]) * maxs;
    double alpha = (-1.0 + 2.0 * params[1]) * maxalpha;

    double z = alpha * s*s;

    double cs[4];
    for(int i = 0; i < 4; ++i)
        cs[i] = universal_var_stumpff_series(i, z);

    for(int i = 0; i < 4; ++i)
        ASSERT(isfinite(cs[i]), "Stumpff series c%d not NaN", i);

    for(int i = 0; i < 2; ++i)
        ASSERT_EQF(1.0 - cs[i], z*cs[i+2], "Stumpff series c%d recurrence relation", i);

    double cs_fun[4]  = {
        universal_var_stumpff_c0(alpha, s),
        universal_var_stumpff_c1(alpha, s),
        universal_var_stumpff_c2(alpha, s),
        universal_var_stumpff_c3(alpha, s)
    };

    for(int i = 0; i < 4; ++i)
        ASSERT(isfinite(cs_fun[i]), "Stumpff function c%d not NaN", i);

    for(int i = 0; i < 2; ++i)
        ASSERT_EQF(1.0 - cs_fun[i], z*cs_fun[i+2], "Stumpff function c%d recurrence relation", i);

    for(int i = 0; i < 4; ++i)
        ASSERT_EQF(cs[i], cs_fun[i], "Stumpff function and series c%d are equal", i);

}

void universal_var_CS_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 1, "");

    const double maxz = 4.0*M_PI;
    double z = (-1.0 + 2.0 * params[0]) * maxz;
    double cf = universal_var_C(z), cs = universal_var_Cseries(z);
    double sf = universal_var_S(z), ss = universal_var_Sseries(z);

    ASSERT(isfinite(cf) && isfinite(cs) && isfinite(sf) && isfinite(ss),
        "C and S not NaN");

    ASSERT_EQF(cs, cf, "C series equals C function w/trigonometry");
    ASSERT_EQF(ss, sf, "S series equals S function w/trigonometry");

    double dCdz = universal_var_dCdz(z, cs, ss);
    double dSdz = universal_var_dSdz(z, cs, ss);

    if(ZEROF(z)) // TODO: derivatives at z = 0
        dCdz = dSdz = 0.0;

    ASSERT(isfinite(dCdz) && isfinite(dSdz),
        "dC/dz and dS/dz not NaN");

    double dz = 1.0e-10 * maxz;
    double Cplus = universal_var_Cseries(z + dz);
    double Cminus = universal_var_Cseries(z - dz);
    double Splus = universal_var_Sseries(z + dz);
    double Sminus = universal_var_Sseries(z - dz);

    ASSERT_EQF(2.0 * dCdz * dz, Cplus - Cminus, "dC/dz");
    ASSERT_EQF(2.0 * dSdz * dz, Splus - Sminus, "dS/dz");
}

void universal_var_from_elements_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 4, "");

    double p = 1.0 + params[0] * 1.0e8;
    double e = params[1] * 5.0;
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

    double k1 = -1.0 + params[3] * 2.0;
    double f1 = maxf * k1;
    double E1 = kepler_anomaly_true_to_eccentric(e, f1);
    double M1 = kepler_anomaly_eccentric_to_mean(e, E1);
    double t1 = t0 + M1 / n;

    double k2 = k1 < 0 ? k1+1.0 : k1-1.0;
    double f2 = maxf * k2;
    double E2 = kepler_anomaly_true_to_eccentric(e, f2);
    double M2 = kepler_anomaly_eccentric_to_mean(e, E2);
    double t2 = t0 + M2 / n;

    double pos1[3], vel1[3], pos2[3], vel2[3];
    kepler_elements_to_state_E(&elements, E1, pos1, vel1);
    kepler_elements_to_state_E(&elements, E2, pos2, vel2);

    double r1 = mag(pos1), v12 = dot(vel1, vel1);
    double rv = dot(pos1, vel1);
    double sqrtmu = sqrt(mu);
    double alpha = (2.0*mu/r1 - v12) / mu;

    double x0 = universal_var_guess_x(mu, sqrtmu, alpha, rv, r1, t2-t1);
    double x = universal_var_iterate_x(sqrtmu, alpha, rv, r1, x0, t2-t1);
    double z = alpha * square(x);

    ASSERT(isfinite(x0),
        "Universal variable x initial guess not NaN");
    ASSERT(isfinite(x),
        "Universal variable x not NaN");
    ASSERT(isfinite(z),
        "Universal variable z not NaN");

    if(kepler_orbit_parabolic(&elements)) {
        ASSERT_EQF(E2-E1, x,
            "Parabolic anomaly and universal variable x = D-D0");
    } else if(kepler_orbit_hyperbolic(&elements)) {
        ASSERT_EQF(-square(E2-E1), z,
            "Hyperbolic anomaly and universal variable z = -(F-F0)^2");
    } else {
        ASSERT_EQF(square(E2-E1), z,
            "Eccentric anomaly and universal variable z = (E-E0)^2");
    }

    double C = universal_var_Cseries(z);
    double S = universal_var_Sseries(z);

    ASSERT(isfinite(C) && isfinite(S),
        "Universal variable C and S not NaN");

    double f, g;
    universal_var_fg(sqrtmu, r1, t2-t1, x, z, C, S, &f, &g);

    ASSERT(isfinite(f) && isfinite(g),
        "Universal variable f and g not NaN");

    double fE, gE;
    universal_var_fg_E(sqrtmu, a, e, r1, t2-t1, E2-E1, &fE, &gE);

    ASSERT(isfinite(fE) && isfinite(gE),
        "Universal variable f(E) and g(E) not NaN");

    ASSERT_EQF(fE, f,
        "Universal var f and eccentric anomaly match");
    ASSERT_EQF(gE, g,
        "Universal var g and eccentric anomaly match");

    double pos3[3];
    for(int i = 0; i < 3; ++i)
        pos3[i] = fE * pos2[i] + gE * vel2[i]; // XXX: use f,g. not fE,gE
    double r3 = mag(pos3);

    double fdot, gdot;
    universal_var_fgdot(sqrtmu, r1, r3, x, z, C, S, &fdot, &gdot);

    ASSERT(isfinite(fdot) && isfinite(gdot),
        "Universal var fdot and gdot not NaN");

    double fdotE, gdotE;
    universal_var_fgdot_E(sqrtmu, a, e, r1, r3, t2-t1, E2-E1, &fdotE, &gdotE);

    ASSERT(isfinite(fdotE) && isfinite(gdotE),
        "Universal var fdot(E) and gdot(E) not NaN");

    ASSERT_EQF(fdot, fdotE,
        "Universal var fdot and eccentric anomaly match");
    ASSERT_EQF(gdot, gdotE,
        "Universal var gdot and eccentric anomaly match");

    double ff, gf;
    universal_var_fg_f(sqrtmu, p, r1, r3, f2-f1, &ff, &gf);
    ASSERT(isfinite(ff) && isfinite(gf),
        "Universal var f(f) and g(f) not NaN");

    double fdotf, gdotf;
    universal_var_fgdot_f(sqrtmu, p, r1, r3, f2-f1, &fdotf, &gdotf);
    ASSERT(isfinite(fdotf) && isfinite(gdotf),
        "Universal var fdot(f) and gdot(f) not NaN");

    ASSERT_EQF(ff, fE,
        "Universal var f true and eccentric anomaly match");
    ASSERT_EQF(gf, gE,
        "Universal var g true and eccentric anomaly match");
    ASSERT_EQF(fdotf, fdotE,
        "Universal var fdot true and eccentric anomaly match");
    ASSERT_EQF(gdotf, gdotE,
        "Universal var g true and eccentric anomaly match");

    ASSERT_EQF(ff, f,
        "Universal var f and true anomaly match");
    ASSERT_EQF(gf, g,
        "Universal var g and true anomaly match");
    ASSERT_EQF(fdotf, fdot,
        "Universal var and fdot true anomaly match");
    ASSERT_EQF(gdotf, gdot,
        "Universal var g and true anomaly match");

    ASSERT_EQF(f*gdot - fdot*g, 1.0,
        "f*gdot - fdot*g = 1.0");
    ASSERT_EQF(fE*gdotE - fdotE*gE, 1.0,
        "f*gdot - fdot*g = 1.0  (E)");
    ASSERT_EQF(ff*gdotf - fdotf*gf, 1.0,
        "f*gdot - fdot*g = 1.0  (f)");

    double pos[3], vel[3];
    universal_var(mu, pos1, vel1, t2-t1, pos, vel);

    ASSERT(isfinite(pos[0]) && isfinite(pos[1]) && isfinite(pos[2]) &&
        isfinite(vel[0]) && isfinite(vel[1]) && isfinite(vel[2]),
        "Position and velocity not NaN");

    ASSERT(eqv3(pos, pos2), "Position equality");
    ASSERT(eqv3(vel, vel2), "Velocity equality");

    double r = mag(pos), v = mag(vel);
    (void)r; (void)v;
}
