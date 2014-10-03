#include <libkepler/kepler.h>
#include <libkepler/intercept.h>

#include <math.h>
#include <float.h>

#include "numtest.h"

static void cross(const double *a, const double *b, double *c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = -(a[0]*b[2] - a[2]*b[0]);
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static double dot(const double *a, const double *b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }

static double mag(const double *a) { return sqrt(dot(a, a)); }

static double sign(double x) { return x < 0 ? -1.0 : 1.0; }
static double square(double x) { return x*x; }
static double cube(double x) { return x*x*x; }

static double clamp(double min, double max, double x) { return (x < min ? min : (x > max ? max : x)); }


static void elliptical_orbit(
    double mu,
    double ap, double pe,
    struct kepler_elements *elements) {
    double a = (ap + pe)/2.0;
    double e = (ap - pe) / (ap + pe);
    double p = a * (1 - e*e);
    double n = sqrt(mu / a*a*a);

    elements->semi_latus_rectum = p;
    elements->eccentricity = e;
    elements->mean_motion = n;
}

void intercept_simple_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
    (void)extra_args;

    ASSERT(num_params == 1, "");

    (void)params;

    // p, e, n, i, an, arg, tau
    double mu = 1.0;
    struct kepler_elements orbit1 = { 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
    struct kepler_elements orbit2 = { 1.0, 0.0, 1.0, M_PI/8.0, 0.0, 0.0, 0.0 };

    elliptical_orbit(mu, 2.0, 1.0, &orbit1);
    elliptical_orbit(mu, 1.0, 1.0, &orbit2);

    struct intercept intercept;
    bool hit = intercept_orbit(&orbit1, &orbit2, -M_PI, M_PI, &intercept);

    ASSERT(hit, "hit");
}

#include <stdio.h>
void print_orbit(const struct kepler_elements *orbit, int idx) {
    printf("p%d=%lf ; e%d=%lf ; i%d=%lf ; an%d=%lf ; argp%d=%lf\n",
        idx, orbit->semi_latus_rectum,
        idx, orbit->eccentricity,
        idx, orbit->inclination,
        idx, orbit->longitude_of_ascending_node,
        idx, orbit->argument_of_periapsis);
}

void intercept_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
    (void)extra_args;

    ASSERT(num_params == 10, "");

    double mu = 1.0 + params[0] * 1.0e10;

    double p1 = 1.0 + params[1] * 1.0e10;
    double e1 = params[2] * 4.0;
    double a1 = p1 / (1.0 - square(e1));
    double n1 = ZEROF(e1 - 1.0) ?
        sqrt(mu / cube(p1)) :    // parabolic
        sqrt(mu / cube(fabs(a1)));

    double i = params[3] * M_PI;
    double an = ZEROF(i) ? 0.0 : (-1.0 + 2.0 * params[4]) * M_PI;
    double arg = ZEROF(e1) ? 0.0 : (-1.0 + 2.0 * params[5]) * M_PI;

    struct kepler_elements orbit1 = { p1, e1, n1, i, an, arg, 0.0 };

    double maxf1 = e1 > 1.0 ?
        acos(1.0/e1) :       // hyperbolic
        (ZEROF(e1-1.0) ?
            4*M_PI/5 :      // parabolic
            M_PI);          // closed orbit
    double f1 = (-1.0 + 2.0 * params[6]) * maxf1;
    double t0 = kepler_orbit_periapsis_time(&orbit1) +
        kepler_anomaly_true_to_mean(e1, f1) / kepler_orbit_mean_motion(&orbit1);

    double pos1[3], vel1[3];
    kepler_elements_to_state_t(&orbit1, t0, pos1, vel1);

    double r1 = mag(pos1);
    double rad1[3] = { pos1[0]/r1, pos1[1]/r1, pos1[2]/r1 };
    double tan1[3], nor1[3];
    kepler_orbit_normal(&orbit1, nor1);
    cross(nor1, rad1, tan1);

    double e2 = params[7] * 4.0;
    double maxf2 = e2 > 1.0 ?
        acos(1.0/e2) :       // hyperbolic
        (ZEROF(e2-1.0) ?
            4*M_PI/5 :      // parabolic
            M_PI);          // closed orbit
    double f2 = (-1.0 + 2.0 * params[8]) * maxf2;

    double p2 = r1 * (1.0 + e2 * cos(f2));
    double a2 = p2 / (1.0 - e2*e2);
    double v2 = ZEROF(e2-1.0) ?
        sqrt(mu * (2.0 / r1)) :
        sqrt(mu * (2.0 / r1 - 1.0 / a2));

    // radial and tangential velocity
    double r_pe2 = p2 / (1 + e2), v_pe2 = sqrt((mu / p2) * square(1.0 + e2));
    double h2 = r_pe2*v_pe2;
    double phi2 = sign(f2) * acos(clamp(-1.0, 1.0, (h2)/(r1*v2)));
    double vt2 = v2 * cos(phi2), vr2 = v2 * sin(phi2);

    double reli = (-1.0 + 2.0 * params[9]) * M_PI;

    double vel2[3] = {
        (cos(reli)*tan1[0] + sin(reli)*nor1[0]) * vt2 + rad1[0]*vr2,
        (cos(reli)*tan1[1] + sin(reli)*nor1[1]) * vt2 + rad1[1]*vr2,
        (cos(reli)*tan1[2] + sin(reli)*nor1[2]) * vt2 + rad1[2]*vr2,
    };

    struct kepler_elements orbit2;
    kepler_elements_from_state(mu, pos1, vel2, t0, &orbit2);

    struct kepler_elements orbits[2] = { orbit1, orbit2 };

    //printf("reli: %lf\n", reli);
    double pos2[3];
    kepler_elements_to_state_t(&orbit2, t0, pos2, vel2);

    // "fix" true anomaly of orbit2 (f2 is wrong for circular orbits)
    double tan2[3], bit2[3];
    kepler_orbit_tangent(&orbit2, tan2);
    kepler_orbit_bitangent(&orbit2, bit2);
    double f22 = sign(dot(bit2, pos1)) * acos(clamp(-1.0, 1.0, dot(pos1, tan2)/mag(pos1)));

    //printf("f2: %lf\tf22: %lf\n", f2, f22);
    //print_orbit(&orbit1, 1);
    //print_orbit(&orbit2, 2);

    for(int i = 0; i < 2; ++i) {
        double fs[4];
        double threshold = (1.0/1000.0) *
            (kepler_orbit_semi_latus_rectum(&orbit1) +
             kepler_orbit_semi_latus_rectum(&orbit2));
        int intersects = intersect_orbit(&orbits[i], &orbits[!i], threshold, fs);

        ASSERT(intersects == 1 || intersects == 2,
            "One or two intersects");

        if(intersects == 0)
            continue; // XXX: this should not happen

        ASSERT(isfinite(fs[0]) && isfinite(fs[1]) &&
            (intersects < 2 || (isfinite(fs[2]) && isfinite(fs[3]))),
            "Intersect true anomaly not NaN");

        for(int j = 0; j < intersects*2; ++j)
            ASSERT_RANGEF(fs[j], -M_PI, M_PI,
                "Intersect true anomaly in range -pi..pi");

        double f = i == 0 ? f1 : f22;

        ASSERT(
            // intersect1
            (fs[0] <= fs[1] && LTF(fs[0], f) && LTF(f, fs[1])) ||
            // intersect 1 near apoapsis
            (fs[0] > fs[1] && (LTF(fs[0], f) || LTF(f, fs[1]))) ||
            (intersects == 2 && (
            // intersect2
                (fs[2] <= fs[3] && LTF(fs[2], f) && LTF(f, fs[3])) ||
            // intersect2 near apoapsis
                (fs[2] > fs[3] && (LTF(fs[2], f) || LTF(f, fs[3]))))),
            "True anomaly within range");

        ASSERT(intersects < 2 || !(fs[0] > fs[1] && fs[2] > fs[3]),
            "True anomaly ranges overlap at apoapsis");
        ASSERT(intersects < 2 || (fs[0] < fs[1]) ||
                (fs[1] < fs[2] && fs[3] < fs[0]),
            "True anomaly ranges overlap");
        ASSERT(intersects < 2 || (fs[2] < fs[3]) ||
                (fs[3] < fs[0] && fs[1] < fs[2]),
            "True anomaly ranges overlap");
        ASSERT(intersects < 2 || !(fs[0] < fs[1] && fs[2] < fs[3]) ||
            fs[1] < fs[2],
            "True anomaly ranges overlap");
    }
}
