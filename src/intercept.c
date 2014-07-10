#include <math.h>
#include <float.h>

#include <libkepler/kepler.h>
#include <libkepler/intercept.h>

static void cross(const double *a, const double *b, double *c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = -(a[0]*b[2] - a[2]*b[0]);
    c[2] = a[0]*b[1] - a[1]*b[0];
}

static double dot(const double *a, const double *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static double mag(const double *a) {
    return sqrt(dot(a, a));
}

static double clamp(double min, double max, double x) {
    return (x < min ? min : (x > max ? max : x));
}

static bool zero(double x) {
    return x*x < DBL_EPSILON;
}

static double sign(double x) { return x < 0 ? -1.0 : 1.0; }
//static double square(double x) { return x*x; }
//static double cube(double x) { return x*x*x; }

static double min(double a, double b) { return a < b ? a : b; }
static double max(double a, double b) { return a > b ? a : b; }

static double angle_clamp(double x0) {
    double x = (x0+M_PI)/(2.0*M_PI);
    return -M_PI + 2.0*M_PI * (x - floor(x));
}

static void swap(double *a, double *b) { double t = *a; *a = *b; *b = t; }

static double true_anomaly_from_radius(const struct kepler_elements *elements, double r) {
    double p = kepler_orbit_semi_latus_rectum(elements);
    double e = kepler_orbit_eccentricity(elements);

    return acos(clamp(-1.0, 1.0, (p/r - 1.0)/e));
}

int intersect_orbit(
    const struct kepler_elements *orbit1,
    const struct kepler_elements *orbit2,
    double threshold,
    double *fs
    ) {
    // Apoapsis-periapsis test early exit
    if((kepler_orbit_closed(orbit1) &&
            kepler_orbit_apoapsis(orbit1) < (kepler_orbit_periapsis(orbit2) - threshold)) ||
        (kepler_orbit_closed(orbit2) &&
             kepler_orbit_apoapsis(orbit2) < (kepler_orbit_periapsis(orbit1) - threshold)))
        return 0;

    // Altitude check
    double f1 = 0.0, f2 = M_PI;
    if(!kepler_orbit_circular(orbit1)) {
        f1 = true_anomaly_from_radius(orbit1, kepler_orbit_periapsis(orbit2) - threshold);

        if(kepler_orbit_closed(orbit2)) {
            f2 = true_anomaly_from_radius(orbit1, kepler_orbit_apoapsis(orbit2) + threshold);
        }
    }

    // Relative inclination and line of nodes
    double nor1[3], nor2[3];
    double nodes[3];
    kepler_orbit_normal(orbit1, nor1);
    kepler_orbit_normal(orbit2, nor2);
    cross(nor2, nor1, nodes); // points to ascending node where orbit1 moves above orbit2
    double N = mag(nodes);
    bool coplanar = zero(N);
    double rel_incl = asin(clamp(-1.0, 1.0, N));

    // Non-coplanar orbits
    double f_an = 0.0, f_dn = M_PI, delta_f = M_PI/2.0;
    if(coplanar) {
        if(zero(f1)) { fs[0] = -f2; fs[1] = f2; return 1; }       // apoapsis
        else if(!(f2 < M_PI)) { fs[0] = f1; fs[2] = -f1; return 1; } // periapsis
        else { fs[0] = -f2; fs[1] = -f1; fs[2] = f1; fs[3] = f2; return 2; }
    } else {
        double tan1[3], bit1[3];
        kepler_orbit_tangent(orbit1, tan1);
        kepler_orbit_bitangent(orbit1, bit1);

        f_an = sign(dot(bit1, nodes)) * acos(clamp(-1.0, 1.0, dot(nodes, tan1)/N));
        f_dn = f_an - sign(f_an) * M_PI;

        double r_pe = kepler_orbit_periapsis(orbit1); // XXX: ap or pe?
        // spherical trigonometry sine law
        delta_f = asin(clamp(-1.0, 1.0, sin(threshold / (2.0*r_pe)) / sin(rel_incl / 2.0)));
    }

    if(zero(f1) && !(f2 < M_PI)) {
        // intersects anywhere on orbit (f = -pi .. pi)
        fs[0] = angle_clamp(min(f_an, f_dn)-delta_f);
        fs[1] = angle_clamp(min(f_an, f_dn)+delta_f);
        fs[2] = angle_clamp(max(f_an, f_dn)-delta_f);
        fs[3] = angle_clamp(max(f_an, f_dn)+delta_f);

        return 2;
    } else if(zero(f1)) {
        // intersect near periapsis (f = -f2 .. f2)
        fs[0] = max(min(f_an, f_dn)-delta_f, -f2);
        fs[1] = min(min(f_an, f_dn)+delta_f, f2);
        fs[2] = max(max(f_an, f_dn)-delta_f, -f2);
        fs[3] = min(max(f_an, f_dn)+delta_f, f2);

        if(fs[1] >= fs[2]) {
            fs[1] = fs[3];
            return 1;
        }
    } else if(!(f2 < M_PI)) {
        // intersect near apoapsis (f < -f1, f > f1)
        fs[0] = angle_clamp(max(min(f_an, f_dn)-delta_f, -f1));
        fs[1] = min(min(f_an, f_dn)+delta_f, -f1);
        fs[2] = max(max(f_an, f_dn)-delta_f, f1);
        fs[3] = angle_clamp(min(max(f_an, f_dn)+delta_f, f1));

        if(fs[0] > 0.0 && fs[3] < 0.0) {
            fs[0] = fs[2];
            return 1;
        }
    } else {
        // two intersects (-f2 < f < -f1, f1 < f < f2)
        fs[0] = angle_clamp(max(min(f_an, f_dn)-delta_f, -f2));
        fs[1] = min(min(f_an, f_dn)+delta_f, -f1);
        fs[2] = max(max(f_an, f_dn)-delta_f, f1);
        fs[3] = angle_clamp(min(max(f_an, f_dn)+delta_f, f2));

        if(fs[0] > 0.0 && fs[3] < 0.0) {
            fs[0] = fs[2];
            return 1;
        }

        if(fs[1] >= fs[2]) {
            fs[1] = fs[3];
            return 1;
        }
    }

    if((fs[0] > fs[1]) && !(fs[2] > fs[3])) { swap(fs+0, fs+2); swap(fs+1, fs+3); }


 /*
    printf("f1: %2.3lf\t f2: %2.3lf\t f_an: %2.3lf\t f_dn: %2.3lf\t delta_f: %2.3lf\n",
        f1, f2, f_an, f_dn, delta_f);
    printf("f11: %2.3lf\t f12: %2.3lf\t f21: %2.3lf\t f22: %2.3lf\n",
        fs[0], fs[1], fs[2], fs[3]);
        */

    return (fs[0] < fs[1]) + (fs[2] < fs[3]);
}

bool intercept_minimize(
    const struct kepler_elements *orbit1,
    const struct kepler_elements *orbit2,
    double threshold,
    double t0, double t1,
    struct intercept *intercept) {
    (void)orbit1;
    (void)orbit2;
    (void)threshold;
    (void)t0;
    (void)t1;
    (void)intercept;

    return false;
}

bool intercept_orbit(
    const struct kepler_elements *orbit1,
    const struct kepler_elements *orbit2,
    double t0, double t1,
    struct intercept *intercept) {
    (void)intercept;

    // TODO: adjustable threshold variable, SOI search
    double threshold = (1.0/1000.0) *
        min(kepler_orbit_semi_latus_rectum(orbit1),
            kepler_orbit_semi_latus_rectum(orbit2));

    const struct kepler_elements *orbits[2] = { orbit1, orbit2 };

    // true anomaly ranges of possible intercepts
    double fs[2][4];
    int is[2];

    for(int i = 0; i < 2; ++i)
        if((is[i] = intersect_orbit(orbits[i], orbits[!i], threshold, &fs[i][0])) == 0)
            return false;

    // time ranges of possible intercepts
    double ts[2][4];
    for(int o = 0; o < 2; ++o)
        for(int i = 0; i < 2*is[o]; ++i)
            ts[o][i] =
                kepler_orbit_periapsis_time(orbits[o]) +
                (kepler_anomaly_true_to_mean(kepler_orbit_eccentricity(orbits[o]), fs[o][i]) /
                 kepler_orbit_mean_motion(orbits[o]));

    double dt = t1 - t0;
    if(kepler_orbit_closed(orbits[0]) || kepler_orbit_closed(orbits[1])) {
        int smaller = (!kepler_orbit_closed(orbits[0]) ||
             kepler_orbit_period(orbits[1]) < kepler_orbit_period(orbits[0]));
        dt = kepler_orbit_period(orbits[smaller]);
    }

    (void)ts;
    for(double t = t0; t < t1; t += dt) {
        //if(kepler_orbit_closed(orbit1))
            //double tp = fmod(t - kepler_orbit_periapsis_time(orbit1), kepler_orbit_period(orbit1);
    }

    return false;
}
