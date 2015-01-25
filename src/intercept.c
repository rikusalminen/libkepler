#include <libkepler/conic.h>
#include <libkepler/anomaly.h>
#include <libkepler/true_anomaly.h>
#include <libkepler/orbit.h>

struct intercept {
    double time;
    vec4d relative_position;
    vec4d relative_velocity;
};

int intercept_intersect(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double threshold,
    double *fs) {
    // radial orbits not handled
    if(orbit_radial(orbit1) || orbit_radial(orbit2))
        return 0;

    // find 0, 1, or 2 ranges of true anomaly (w.r.t. orbit1)
    // where orbit1 and orbit2 may be than threshold

    double p1 = orbit_semi_latus_rectum(orbit1);
    double p2 = orbit_semi_latus_rectum(orbit2);
    double e1 = orbit_eccentricity(orbit1);
    double e2 = orbit_eccentricity(orbit2);

    // apoapsis-periapsis test
    if((conic_closed(e1) && conic_apoapsis(p1, e1) < conic_periapsis(p2, e2)) ||
        (conic_closed(e2) && conic_apoapsis(p2, e2) < conic_periapsis(p1, e1)))
        return 0;

    // true anomaly between target apoapsis/periapsis +/- threshold
    double maxf = conic_hyperbolic(e1) ? acos(1.0 / e1) : M_PI;
    double f1 = conic_circular(e1) ? 0.0 :
        true_anomaly_from_radius(p1, e1, conic_periapsis(p2, e2) - threshold);
    double f2 = fmin(maxf, (conic_circular(e1) || !conic_closed(e2)) ? M_PI :
        true_anomaly_from_radius(p1, e1, conic_apoapsis(p2, e2) + threshold));

    // ascending node
    vec4d nodes = cross(orbit1->normal_axis, orbit2->normal_axis);
    double N2 = dot(nodes, nodes)[0], N = sqrt(N2);
    int coplanar = N2 < DBL_EPSILON;
    double reli = sign(dot(orbit1->normal_axis, orbit2->normal_axis)[0]) *
        asin(clamp(-1.0, 1.0, N));

     // coplanar orbits
    if(coplanar) {
        if(zero(f1)) {
            // intersect near periapsis
            fs[0] = -f2; fs[1] = f2;
            return 1;
        } else if(conic_closed(e1) && !(f2 < M_PI)) {
            // intersect near apoapsis
            fs[0] = f1; fs[1] = -f1;
            return 1;
        } else {
            fs[0] = -f2; fs[1] = -f1; fs[2] = f1; fs[3] = f2;
            return 2;
        }
    }

    // non-coplanar orbits
    double f_an = sign(dot(orbit1->minor_axis, nodes)[0]) *
        acos(clamp(-1.0, 1.0, dot(orbit1->major_axis, nodes)/N));
    double f_dn = f_an - sign(f_an) * M_PI;

    double f_nodes[2] = { fmin(f_an, f_dn), fmax(f_an, f_dn) };

    // TODO: handle non-coplanar orbits
}

int intercept_times(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double t0, double t1,
    const double *fs
    double *times,
    int max_times) {
    // find at most max_times time ranges between (t0..t1)
    // where orbit1 and orbit2 are between (fs[0]..fs[1]) and (fs[2]..fs[3])
}

int intercept_search(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double t0, double t1,
    double threshold,
    double *times) {
    // find 0, 1, or 2 time ranges
    // where distance between orbit1 and orbit2 is less than threshold
}

double intercept_minimize(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double t0, double t1,
    double target_distance,
    struct intercept *intercept) {
    // find 1 time
    // when distance between orbit1 and orbit2 is closest to target_distance
    // there must be an intercept between t0 and t1
}

int intercept_orbit(
    const struct orbit *orbit1,
    const struct orbit *orbit2,
    double t0, double t1,
    double threshold,
    double target_distance,
    struct intercept *intercepts,
    int max_intercepts) {
    // find at most max_intercepts times
    // when distance between orbit1 and orbit2 is closest to target_distance

    // geometry prefilter
    double fs[4];
    int intersects = intercept_intersect(orbit1, orbit2, threshold, fs);
    if(intersects == 0)
        return 0;

    int num_intercepts = 0;
    while(t0 < t1 && num_intercepts < max_intercepts) {
        // time prefilter
        int max_times = max_intercepts;
        double times[2*max_times];
        int num_times = intercept_times(
            orbit1, orbit2,
            t0, t1,
            fs,
            times, max_times);

        if(num_times == 0)
            break;

        for(int i = 0; i < num_times; ++i) {
            double ti[4];
            int is = intercept_search(
                orbit1, orbit2,
                times[2*i+0], times[2*i+1],
                threshold, ti);

            for(int j = 0; j < is; ++i) {
                intercept_minimize(
                    orbit1, orbit2,
                    ti[2*j+0], ti[2*j+1],
                    target_distance,
                    intercepts + num_intercepts);
                t0 = ti[2*j+1];
            }

            num_intercepts += is;
        }
    }

    return num_intercepts;
}
