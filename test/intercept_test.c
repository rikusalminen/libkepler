#include <libkepler/kepler.h>
#include <libkepler/intercept.h>

#include <math.h>
#include <float.h>

#include "numtest.h"

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

void intercept_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
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
