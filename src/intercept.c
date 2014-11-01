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
static double square(double x) { return x*x; }
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

    if(kepler_orbit_hyperbolic(orbit1)) {
        f2 = acos(1.0 / kepler_orbit_eccentricity(orbit1));
    }

    // Relative inclination and line of nodes
    double nor1[3], nor2[3];
    double nodes[3];
    kepler_orbit_normal(orbit1, nor1);
    kepler_orbit_normal(orbit2, nor2);
    cross(nor2, nor1, nodes); // points to ascending node where orbit1 moves above orbit2
    double N = mag(nodes);
    bool coplanar = zero(N);
    double rel_incl = sign(dot(nor1, nor2)) * asin(clamp(-1.0, 1.0, N));

    // Coplanar orbits
    if(coplanar) {
        if(zero(f1)) {
            // intersect near periapsis
            fs[0] = -f2; fs[1] = f2;
            return 1;
        } else if(kepler_orbit_closed(orbit1) && !(f2 < M_PI)) {
            // intersect near apoapsis
            fs[0] = f1; fs[1] = -f1;
            return 1;
        } else {
            fs[0] = -f2; fs[1] = -f1; fs[2] = f1; fs[3] = f2;
            return 2;
        }
    }

    // Non-coplanar orbits
    double tan1[3], bit1[3];
    kepler_orbit_tangent(orbit1, tan1);
    kepler_orbit_bitangent(orbit1, bit1);

    double f_an = sign(dot(bit1, nodes)) * acos(clamp(-1.0, 1.0, dot(nodes, tan1)/N));
    double f_dn = f_an - sign(f_an) * M_PI;

    double f_nodes[2] = { min(f_an, f_dn), max(f_an, f_dn) };
    double delta_fs[2] = { M_PI, M_PI };

    for(int i = 0; i < 2; ++i) {
        // distance at node
        //double r = kepler_orbit_semi_latus_rectum(orbit1) /
            //(1.0 + kepler_orbit_eccentricity(orbit1) * cos(f_nodes[i]));
        double r = kepler_orbit_periapsis(orbit1); // XXX: periapsis?
        // spherical trigonometry sine law
        delta_fs[i] = asin(clamp(-1.0, 1.0, sin(threshold / (2.0*r)) / sin(fabs(rel_incl) / 2.0)));
    }

    if(kepler_orbit_closed(orbit1) && zero(f1) && !(f2 < M_PI)) {
        // intersects anywhere on orbit (f = -pi .. pi)
        fs[0] = angle_clamp(f_nodes[0] - delta_fs[0]);
        fs[1] = angle_clamp(f_nodes[0] + delta_fs[0]);
        fs[2] = angle_clamp(f_nodes[1] - delta_fs[1]);
        fs[3] = angle_clamp(f_nodes[1] + delta_fs[1]);

        return 2;
    } else if(zero(f1)) {
        // intersect near periapsis (f = -f2 .. f2)
        fs[0] = max(f_nodes[0] - delta_fs[0], -f2);
        fs[1] = min(f_nodes[0] + delta_fs[0], f2);
        fs[2] = max(f_nodes[1] - delta_fs[1], -f2);
        fs[3] = min(f_nodes[1] + delta_fs[1], f2);

        if(fs[1] >= fs[2]) {
            fs[1] = fs[3];
            return 1;
        }
    } else if(kepler_orbit_closed(orbit1) && !(f2 < M_PI)) {
        // intersect near apoapsis (f < -f1, f > f1)
        int is = 0;
        if(f_nodes[0] - delta_fs[0] < -f1) {
            fs[2*is+0] = f_nodes[0] - delta_fs[0] < -M_PI ?
                max(angle_clamp(f_nodes[0]-delta_fs[0]), f1) :
                f_nodes[0] - delta_fs[0];
            fs[2*is+1] = min(f_nodes[0] + delta_fs[0], -f1);
            is += 1;
        }

        if(f_nodes[1] + delta_fs[1] > f1) {
            fs[2*is+0] = max(f_nodes[1] - delta_fs[1], f1);
            fs[2*is+1] = f_nodes[1] + delta_fs[1] > M_PI ?
                min(angle_clamp(f_nodes[1] + delta_fs[1]), -f1) :
                f_nodes[1] + delta_fs[1];
            is += 1;
        }

        return is;
    } else {
        // two intersects (-f2 < f < -f1, f1 < f < f2)
        fs[0] = max(f_nodes[0] - delta_fs[0], -f2);
        fs[1] = min(f_nodes[0] + delta_fs[0], -f1);
        fs[2] = max(f_nodes[1] - delta_fs[1], f1);
        fs[3] = min(f_nodes[1] + delta_fs[1], f2);

        if(fs[1] >= fs[2]) {
            fs[1] = fs[3];
            return 1;
        }
    }

    if((fs[0] > fs[1]) && !(fs[2] > fs[3])) {
        swap(fs+0, fs+2);
        swap(fs+1, fs+3);
    }

    return (fs[0] < fs[1]) + (fs[2] < fs[3]);
}

bool intercept_minimize(
    const struct kepler_elements *orbit1,
    const struct kepler_elements *orbit2,
    double threshold,
    double t0, double t1,
    struct intercept *intercept) {
    (void)intercept;

    // this algorithm finds 0, 1 or 2 intercepts during interval t0..t1

    double mu = kepler_orbit_gravity_parameter(orbit1);

    double vmax = kepler_orbit_periapsis_vel(orbit1) +
        kepler_orbit_periapsis_vel(orbit2);

    double amax = mu / square(kepler_orbit_periapsis(orbit1)) +
        mu / square(kepler_orbit_periapsis(orbit2));

    double times[4] = { t0, t1, 0.0, 0.0 };
    int intercepts = 0;

    for(int isect = 0; isect == 0; ++isect) { // XXX: condition

        int converged = 0;
        for(int i = 0; times[2*isect+1] > times[2*isect+0] && converged != 0x3; i = !i) {
            if(converged & (1 << i))
                continue;

            double t = times[2*isect + i];
            double neg = i ? 1.0 : -1.0;

            double pos1[3], vel1[3], pos2[3], vel2[3];
            kepler_elements_to_state_t(orbit1, t, pos1, vel1);
            kepler_elements_to_state_t(orbit2, t, pos2, vel2);

            double dr[3], dv[3];
            for(int j = 0; j < 3; ++j) {
                dr[j] = pos2[j] - pos1[j];
                dv[j] = vel1[j] - vel2[j];
            }

            double d2 = dot(dr, dr), d = sqrt(d2);

            if(fabs(d) < threshold) {
                converged |= 1 << i;
            }

            double v_rel = dot(dr, dv) / d;
            double dt_vel = v_rel < 0 ? -1.0 : d/v_rel;

            if(v_rel * neg > 0.0) {
                // getting closer
            } else {
                // going away
            }

            double a = 0.5 * amax;
            double b = vmax;
            double c = -(fabs(d) - threshold);
            double dt_max = (-b + sqrt(b*b - 4.0*a*c))/(2.0 * a);

            double dt = max(0.0, max(dt_max, dt_vel));
            times[2*isect + i] = t + neg * dt;
        }

        if(converged != 0x3) {
            return false;
        }

        double threshold_t = threshold / vmax;
        if(times[2*isect+1] - times[2*isect+0] > threshold_t) {
            // XXX: two intercepts, split range in two
            //times[2*isect+1]
        }
    }

    return intercepts != 0; // XXX: return value
}

bool intercept_orbit(
    const struct kepler_elements *orbit1,
    const struct kepler_elements *orbit2,
    double t0, double t1,
    struct intercept *intercept) {
    // TODO: adjustable threshold variable, SOI search
    double threshold = (1.0/1000.0) *
        min(kepler_orbit_semi_latus_rectum(orbit1),
            kepler_orbit_semi_latus_rectum(orbit2));

    const struct kepler_elements *orbits[2] = { orbit1, orbit2 };

    // true anomaly ranges of possible intercepts
    double fs[2][4];
    int intersects[2];

    for(int o = 0; o < 2; ++o)
        if((intersects[o] = intersect_orbit(orbits[o], orbits[!o], threshold, &fs[o][0])) == 0)
            return false;

    // time ranges of possible intercepts
    double times[2][4];
    for(int o = 0; o < 2; ++o) {
        double e = kepler_orbit_eccentricity(orbits[o]);
        double n = kepler_orbit_mean_motion(orbits[o]);
        double P = kepler_orbit_period(orbits[o]);
        double t_pe = kepler_orbit_periapsis_time(orbits[o]);

        // intersect #2 is before #1 (near apoapsis)
        int swapped = intersects[o] == 2 && (fs[o][2] > fs[o][3]);

        for(int i = 0; i < intersects[o]; ++i) {
            int isect = swapped ? !i : i;

            for(int j = 0; j < 2; ++j) {
                double f = fs[o][2*isect+j];
                double M = kepler_anomaly_true_to_mean(e, f);

                times[o][2*i+j] = t_pe + (M / n);
            }

            if(times[o][2*i] > times[o][2*i+1]) // intersect near apoapsis
                times[o][2*i] -= P;
        }
    }

    // number of orbital periods to periapsis
    int n_orbit[2] = { 0, 0 };
    int isect[2] = { 0, 0 };

    for(int o = 0; o < 2; ++o) {
        if(kepler_orbit_closed(orbits[o])) {
            double t_pe = kepler_orbit_periapsis_time(orbits[o]);
            double P = kepler_orbit_period(orbits[o]);

            n_orbit[o] = (int)trunc((t0 - t_pe)/P + (t_pe > t0 ? -0.5 : 0.5));
        }
    }

    // loop over orbits to cover entire time range
    int num_intercepts = 0;
    double t = t0;
    while(t < t1) {
        double trange[2][2];

        // time interval on this period
        for(int o = 0; o < 2; ++o) {
            double period = (!kepler_orbit_closed(orbits[o]) ? 0 :
                 n_orbit[o] * kepler_orbit_period(orbits[o]));

            trange[o][0] = times[o][2*isect[o]] + period;
            trange[o][1] = times[o][2*isect[o]+1] + period;
        }

        // overalap of two time intervals
        double t_begin = max(t, max(trange[0][0], trange[1][0]));
        double t_end = min(t1, min(trange[0][1], trange[1][1]));
        t = t_end;

        if(t_begin < t_end) {
            // possible intercept on time interval t_begin .. t_end

            intercept_minimize(orbit1, orbit2, threshold, t_begin, t_end, intercept);

            num_intercepts += 1; // XXX:
        }

        // advance to next time range on one orbit
        int advance = trange[0][1] < trange[1][1] ? 0 : 1;
        isect[advance] += 1;

        if(isect[advance] == intersects[advance]) {
            if(!kepler_orbit_closed(orbits[advance]))
                break;

            // advance to next orbital period
            isect[advance] = 0;
            n_orbit[advance] += 1;
        }
    }

    return num_intercepts != 0; // XXX: return value!
}
