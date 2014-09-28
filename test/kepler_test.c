#include <math.h>
#include <stdbool.h>

#include <libkepler/kepler.h>

#include "numtest.h"

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

static double sign(double x) { return x < 0 ? -1.0 : 1.0; }
static double square(double x) { return x*x; }
static double cube(double x) { return x*x*x; }

static double clamp(double min, double max, double x) {
    return (x < min ? min : (x > max ? max : x));
}

static void matrix_vector_product(const double *m, const double *v, double *x) {
    for(int i = 0; i < 3; ++i) {
        x[i] = 0.0;

        for(int j = 0; j < 3; ++j) {
            x[i] += m[3 * i + j] * v[j];
            //x[i] += m[3 * j + i] * v[j];
        }
    }
}


static bool eqv3(double *a, double *b) {
    double diff[3] = { a[0]-b[0], a[1]-b[1], a[2]-b[2] };
    return (ZEROF(dot(a,a)) && ZEROF(dot(b,b))) ||
        ZEROF(dot(diff, diff)/(dot(a, a) + dot(b, b)));
}

static void rotation_matrix_euler(double rx, double ry, double rz, double *m)
{
    double sn[3] = { sin(rx), sin(ry), sin(rz) };
    double cs[3] = { cos(rx), cos(ry), cos(rz) };

    m[0] = cs[1]*cs[2];
    m[1] = cs[1]*sn[2];
    m[2] = -sn[1];

    m[3] = sn[0]*sn[1]*cs[2] - cs[0]*sn[2];
    m[4] = sn[0]*sn[1]*sn[2] + cs[0]*cs[2];
    m[5] = sn[0]*cs[1];

    m[6] = cs[0]*sn[1]*cs[2] + sn[0]*sn[2];
    m[7] = cs[0]*sn[1]*sn[2] - sn[0]*cs[2];
    m[8] = cs[0]*cs[1];
}

void anomaly_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
    (void)extra_args;

    ASSERT(num_params == 2, "num_params");

    double e = params[0] * 4.0;
    double t = -1.0 + params[1] * 2.0;

    double maxf = e > 1.0 ?
        acos(1.0/e) :       // hyperbolic
        (ZEROF(e-1.0) ?
            4*M_PI/5 :      // parabolic
            M_PI);          // closed orbit
    double maxE = kepler_anomaly_true_to_eccentric(e, maxf);
    double maxM = kepler_anomaly_eccentric_to_mean(e, maxE);

    double M = t * maxM;
    double E = t * maxE;
    double f = t * maxf;

    double f2 = kepler_anomaly_eccentric_to_true(e, E);
    ASSERT_RANGEF(f2, -maxf, maxf, "True anomaly within range");
    ASSERT_EQF(E, kepler_anomaly_true_to_eccentric(e, f2), "Eccentric -> True");

    double E2 = kepler_anomaly_true_to_eccentric(e, f);
    ASSERT_RANGEF(E2, -maxE, maxE, "Eccentric anomaly within range");
    ASSERT_EQF(f, kepler_anomaly_eccentric_to_true(e, E2), "True -> Eccentric");

    double E3 = kepler_anomaly_mean_to_eccentric(e, M);
    ASSERT_RANGEF(E3, -maxE, maxE, "Eccentric anomaly within range");
    ASSERT_EQF(M, kepler_anomaly_eccentric_to_mean(e, E3), "True -> Eccentric");

    double M2 = kepler_anomaly_eccentric_to_mean(e, E);
    ASSERT_RANGEF(M2, -maxM, maxM, "Mean anomaly within range");
    ASSERT_EQF(E, kepler_anomaly_mean_to_eccentric(e, M2), "Eccentric -> Mean");

    double M3 = kepler_anomaly_true_to_mean(e, f);
    ASSERT_RANGEF(M3, -maxM, maxM, "Mean anomaly within range");
    ASSERT_EQF(f, kepler_anomaly_mean_to_true(e, M3), "True -> Mean");

    double f3 = kepler_anomaly_mean_to_true(e, M);
    ASSERT_RANGEF(f3, -maxf, maxf, "True anomaly within range");
    ASSERT_EQF(M, kepler_anomaly_true_to_mean(e, f3), "Mean -> True");

    double dEdM = kepler_anomaly_dEdM(e, E);
    double MM = kepler_anomaly_eccentric_to_mean(e, E);
    double dM = 1.0e-10 * maxM;
    double Eplus = kepler_anomaly_mean_to_eccentric(e, MM+dM);
    double Eminus = kepler_anomaly_mean_to_eccentric(e, MM-dM);
    ASSERT_EQF(dEdM * 2.0 * dM, (Eplus-Eminus), "dE/dM");
}

void orbit_from_state_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
    (void)extra_args;

    ASSERT(num_params == 7, "");

    struct kepler_elements elements = { 0, 0, 0, 0, 0, 0, 0};

    double mu = 1.0 + params[0] * 1.0e10;
    double r0 = 1.0 + params[1] * 1.0e10;
    double v_circ = sqrt(mu / r0), v_esc = sqrt(2.0 * mu / r0);
    double v0 = v_circ + (v_esc - v_circ) * 2.0 * params[2];

    double h0 = r0 * v0;
    double visviva0 = v0*v0/2.0 - mu/r0;
    double e2 = 1.0 + (2.0 * visviva0 * h0*h0)/(mu*mu);
    double e = ZEROF(e2) ? 0.0 : sqrt(e2);
    double p = r0 * (1.0 + e);

    double mat[9];
    double rot[3];
    for(int i = 0; i < 3; ++i)
        rot[i] = -M_PI + 2.0*M_PI * params[3+i];
    rotation_matrix_euler(rot[0], rot[1], rot[2], mat);

    double maxf = e > 1.0 ?
        acos(1.0/e) :       // hyperbolic
        (ZEROF(e-1.0) ?
            4*M_PI/5 :      // parabolic
            M_PI);          // closed orbit
    double f0 = maxf * (-1.0 + 2.0 * params[6]);

    double r = p / (1.0 + e * cos(f0));
    double a = p / (1.0 - e2);
    double v = ZEROF(e-1.0) ?
        sqrt(mu * (2.0 / r)) :
        sqrt(mu * (2.0 / r - 1.0 / a));

    double rx = r * cos(f0), ry = r * sin(f0);

    // flight path angle
    double phi = sign(f0) * acos(clamp(-1.0, 1.0, (r0*v0)/(r*v)));
    double vt = v * cos(phi), vr = v * sin(phi);
    double vx = vr * cos(f0) + vt * -sin(f0);
    double vy = vr * sin(f0) + vt * cos(f0);

    double pos[3] = {
        rx * mat[0] + ry * mat[3],
        rx * mat[1] + ry * mat[4],
        rx * mat[2] + ry * mat[5]
    };

    double vel[3] = {
        vx * mat[0] + vy * mat[3],
        vx * mat[1] + vy * mat[4],
        vx * mat[2] + vy * mat[5]
    };

    double t0 = 0.0;

    kepler_elements_from_state(mu, pos, vel, t0, &elements);

    ASSERT(
        isfinite(kepler_orbit_semi_latus_rectum(&elements)) &&
        isfinite(kepler_orbit_eccentricity(&elements)) &&
        isfinite(kepler_orbit_mean_motion(&elements)) &&
        isfinite(kepler_orbit_inclination(&elements)) &&
        isfinite(kepler_orbit_longitude_of_ascending_node(&elements)) &&
        isfinite(kepler_orbit_argument_of_periapsis(&elements)) &&
        isfinite(kepler_orbit_periapsis_time(&elements)),
        "All elements not NaN or Inf");

    ASSERT_EQF(kepler_orbit_gravity_parameter(&elements), mu,
        "Gravity parameter matches");

    ASSERT_RANGEF(elements.inclination, 0.0, M_PI,
        "Inclination within range 0..pi");
    ASSERT_RANGEF(elements.longitude_of_ascending_node, -M_PI, M_PI,
        "Longitude of ascending node with range -pi..pi");
    ASSERT_RANGEF(elements.argument_of_periapsis, -M_PI, M_PI,
        "Argument of periapsis with range -pi..pi");

    ASSERT_EQF(kepler_orbit_semi_latus_rectum(&elements), p,
        "Semi-latus rectum matches");
    ASSERT((ZEROF(e2) && ZEROF(pow(kepler_orbit_eccentricity(&elements), 2))) ||
        EQF(kepler_orbit_eccentricity(&elements), e),
        "Eccentricity matches");
    ASSERT_EQF(kepler_orbit_specific_angular_momentum(&elements), h0,
        "Specific relative angular momentum matches");
    if(!kepler_orbit_parabolic(&elements)) // XXX: rounding errors near escape
        ASSERT_EQF(kepler_orbit_specific_orbital_energy(&elements), visviva0,
            "Specific orbital energy matches");

    ASSERT(elements.eccentricity >= 0.0,
        "Eccentricity is greater than zero");

    ASSERT(v0 >= v_esc || kepler_orbit_closed(&elements),
        "Closed orbit if less than escape velocity");

    ASSERT_LTF(kepler_orbit_periapsis(&elements), r,
        "Distance is greater than periapsis");
    ASSERT_LTF(v, kepler_orbit_periapsis_vel(&elements),
        "Velocity is smaller than periapsis velocity");

    ASSERT(kepler_orbit_closed(&elements) ==
        !(kepler_orbit_parabolic(&elements) || kepler_orbit_hyperbolic(&elements)),
        "Closed orbits are not parabolic or hyperbolic");
    ASSERT(kepler_orbit_closed(&elements) ||
        kepler_orbit_parabolic(&elements) != kepler_orbit_hyperbolic(&elements),
        "Parabolic orbits are not hyperbolic");

    if(!kepler_orbit_parabolic(&elements)) {
        double visviva = v*v/2.0 - mu/r;
        ASSERT_EQF(visviva, kepler_orbit_specific_orbital_energy(&elements),
            "Specific orbital energy (vis-viva)");
    } else {
        ASSERT(ZEROF(kepler_orbit_specific_orbital_energy(&elements)),
            "Parabolic orbit specific orbital energy is zero (vis-viva)");
    }

    double h[3];
    cross(pos, vel, h);
    ASSERT_EQF(mag(h), kepler_orbit_specific_angular_momentum(&elements),
        "Specific relative angular momentum");

    double pos_pe[3], vel_pe[3];
    kepler_elements_to_state_t(&elements, kepler_orbit_periapsis_time(&elements), pos_pe, vel_pe);

    ASSERT_EQF(mag(pos_pe), kepler_orbit_periapsis(&elements),
        "Periapsis radius");
    ASSERT_EQF(mag(vel_pe), kepler_orbit_periapsis_vel(&elements),
        "Periapsis velocity");

    if(kepler_orbit_closed(&elements)) {
        ASSERT_LTF(kepler_orbit_periapsis(&elements), kepler_orbit_apoapsis(&elements),
            "Apoapsis is greater than periapsis");
        ASSERT_LTF(kepler_orbit_apoapsis_vel(&elements), kepler_orbit_periapsis_vel(&elements),
            "Periapsis velocity is greater than apoapsis velocity");

        ASSERT(kepler_orbit_semi_major_axis(&elements) > 0 && kepler_orbit_semi_minor_axis(&elements) > 0,
            "Closed orbit semi-major and semi-minor axes are positive");
        ASSERT_LTF(kepler_orbit_semi_minor_axis(&elements), kepler_orbit_semi_major_axis(&elements),
            "Semi-minor axis is less than semi_major axis");

        ASSERT(isfinite(kepler_orbit_period(&elements)),
            "Closed orbits have finite period");

        ASSERT_RANGEF(r, kepler_orbit_periapsis(&elements), kepler_orbit_apoapsis(&elements),
            "Distance is between apoapsis and periapsis");
        ASSERT_RANGEF(v, kepler_orbit_apoapsis_vel(&elements), kepler_orbit_periapsis_vel(&elements),
            "Velocity is between apoapsis and periapsis velocity");

        ASSERT(0.0 > kepler_orbit_specific_orbital_energy(&elements),
            "Closed orbits have negative specific orbital energy");

        double pos_ap[3], vel_ap[3];
        kepler_elements_to_state_t(&elements,
                kepler_orbit_periapsis_time(&elements) + kepler_orbit_period(&elements)/2.0,
                pos_ap, vel_ap);

        ASSERT_EQF(mag(pos_ap), kepler_orbit_apoapsis(&elements),
            "Apoapsis radius");
        ASSERT_EQF(mag(vel_ap), kepler_orbit_apoapsis_vel(&elements),
            "Apoapsis velocity");
    }

    if(kepler_orbit_circular(&elements)) {
        ASSERT_EQF(kepler_orbit_apoapsis(&elements), kepler_orbit_periapsis(&elements),
            "Circular orbit periapsis == apoapsis");
        ASSERT_EQF(kepler_orbit_apoapsis_vel(&elements), kepler_orbit_periapsis_vel(&elements),
            "Circular orbit periapsis velocity == apoapsis velocity");
    }

    if(kepler_orbit_parabolic(&elements)) {
        ASSERT(ZEROF(kepler_orbit_specific_orbital_energy(&elements)),
            "Parabolic orbits have zero orbital energy");

        ASSERT_EQF(2.0 * kepler_orbit_periapsis(&elements), kepler_orbit_semi_latus_rectum(&elements),
            "Parabolic orbit periapsis is half semi latus rectum");
    }

    if(kepler_orbit_hyperbolic(&elements)) {
        ASSERT(kepler_orbit_semi_major_axis(&elements) < 0 && kepler_orbit_semi_minor_axis(&elements) < 0,
            "Hyperbolic orbit semi-major and semi-minor axes are positive");

        ASSERT(0.0 < kepler_orbit_specific_orbital_energy(&elements),
            "Hyperbolic orbits have positive specific orbital energy");
    }

    double t[3], n[3], b[3];
    kepler_orbit_tangent(&elements, t);
    kepler_orbit_normal(&elements, n);
    kepler_orbit_bitangent(&elements, b);

    ASSERT(ZEROF(dot(t, n)) && ZEROF(dot(t, b)) && ZEROF(dot(n, b)),
        "Orbit normal, tangent, binormal are orthogonal");

    ASSERT(ZEROF(dot(pos, n)/r) && ZEROF(dot(vel, n)/v),
        "Orbit motion is planar");

    ASSERT_EQF(dot(n, h), mag(h),
        "Orbit normal and angular momentum");

    double matrix[9];
    kepler_orbit_matrix(&elements, matrix);

    for(int pass = 0; pass < 2; ++pass)
    {
        double pos2[3], vel2[3];
        double M = kepler_orbit_mean_anomaly_at_time(&elements, t0);
        double E = kepler_anomaly_mean_to_eccentric(elements.eccentricity, M);
        double f = kepler_anomaly_eccentric_to_true(elements.eccentricity, E);

        ASSERT(ZEROF(e2) ||
            (ZEROF(f0*f0) && ZEROF(f*f)) ||
            (ZEROF(fabs(f0)-M_PI) && ZEROF(fabs(f)-M_PI)) ||
            EQF(f, f0),
            "True anomaly identity");

        if(pass == 0)
            kepler_elements_to_state_E(&elements, E, pos2, vel2);
        else
            kepler_elements_to_state_f(&elements, f, pos2, vel2);

        double pos3[3], vel3[3];
        matrix_vector_product(matrix, pos2, pos3);
        matrix_vector_product(matrix, vel2, vel3);

        ASSERT(eqv3(pos, pos3), "Position identity");
        ASSERT(eqv3(vel, vel3), "Velocity identity");
    }
}

void orbit_from_elements_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
    (void)extra_args;

    ASSERT(num_params == 4, "");

    double p = 1.0 + params[0] * 1.0e10;
    double e = params[1] * 4.0;
    double mu = 1.0 + params[2] * 1.0e10;
    double a = p / (1.0 - square(e));
    double n = ZEROF(e - 1.0) ?
        sqrt(mu / cube(p)) :    // parabolic
        sqrt(mu / cube(fabs(a)));
    double t0 = 0.0;

    struct kepler_elements elements = { p, e, n, 0.0, 0.0, 0.0, t0 };

    ASSERT_EQF(mu, kepler_orbit_gravity_parameter(&elements),
        "Gravitational constant is correct");

    double maxf = e > 1.0 ?
        acos(1.0/e) :       // hyperbolic
        (ZEROF(e-1.0) ?
            4*M_PI/5 :      // parabolic
            M_PI);          // closed orbit

    double x1 = -1.0 + params[3] * 2.0;
    double f1 = maxf * x1;
    double E1 = kepler_anomaly_true_to_eccentric(e, f1);
    double t1 = kepler_orbit_periapsis_time(&elements) +
        kepler_anomaly_eccentric_to_mean(e, E1) / n;

    const int states_size = 2 * 2 * 3; // 2 x position,velocity x vec3
    double state_data[states_size];

    kepler_elements_to_state_f(&elements, f1, &state_data[6*0 + 0], &state_data[6*0 + 3]);
    kepler_elements_to_state_E(&elements, E1, &state_data[6*1 + 0], &state_data[6*1 + 3]);

    for(double *ptr = state_data+0; ptr != state_data + states_size; ++ptr)
        ASSERT(isfinite(*ptr), "Position and velocity not NaN");

    ASSERT(eqv3(state_data + 6*0 + 0, state_data + 6*1 + 0),
        "Conic and eccentric positions match");

    ASSERT(eqv3(state_data + 6*0 + 3, state_data + 6*1 + 3),
        "Conic and eccentric velocities match");

    for(int pass = 0; pass < 2; ++pass) {
        double *ptr = state_data + pass*2*3;
        double *pos = ptr + 0, *vel = ptr + 3;
        double r = mag(pos), v2 = dot(vel, vel); //, v = sqrt(v2);

        ASSERT(dot(pos, pos) > 0 && dot(vel, vel) > 0,
            "Position and Velocity non zero");

        double h[3];
        cross(pos, vel, h);

        ASSERT_EQF(mag(h), kepler_orbit_specific_angular_momentum(&elements),
            "Specific relative angular momentum");

        if(!kepler_orbit_parabolic(&elements)) {
            double visviva = v2/2 - mu/r;
            ASSERT_EQF(visviva, kepler_orbit_specific_orbital_energy(&elements),
                "Specific orbital energy (vis-viva)");
        } else {
            ASSERT_EQF(v2/2, mu/r,
                "Specific orbital energy (vis-viva)");
        }

        if(kepler_orbit_parabolic(&elements)) {
            ASSERT_EQF(
                pos[0] - p/2.0,
                -1.0/(2.0 * p) * square(pos[1]),
                "Orbit is a parabola");

            ASSERT_EQF(dot(pos, vel) / sqrt(mu * p), E1,
                "Parabolic anomaly");
        } else if(kepler_orbit_hyperbolic(&elements)) {
            ASSERT_EQF(1.0,
                square(pos[0] + a*e)/square(kepler_orbit_semi_major_axis(&elements)) -
                square(pos[1])/square(kepler_orbit_semi_minor_axis(&elements)),
                "Orbit is a hyperbola");
        } else {
            ASSERT_EQF(1.0,
                square(pos[0] + a*e)/square(kepler_orbit_semi_major_axis(&elements)) +
                square(pos[1])/square(kepler_orbit_semi_minor_axis(&elements)),
                "Orbit is an ellipse");
        }

        double epoch = elements.periapsis_time +
            kepler_anomaly_true_to_mean(e, f1) / n;

        struct kepler_elements e2;
        kepler_elements_from_state(mu, pos, vel, epoch, &e2);
        ASSERT_EQF(elements.semi_latus_rectum, e2.semi_latus_rectum,
            "Semi-latus rectum identity");
        ASSERT_EQF(elements.eccentricity, e2.eccentricity,
            "Eccentricity identity");
        ASSERT_EQF(elements.mean_motion, e2.mean_motion,
            "Mean motion identity");

        ASSERT(ZEROF(e2.inclination) &&
            ZEROF(e2.longitude_of_ascending_node) &&
            ZEROF(e2.argument_of_periapsis),
            "Zero inclination, longitude of ascending node and argument of periapsis");

        if(kepler_orbit_closed(&elements)) {
            ASSERT(kepler_orbit_closed(&e2),
                "Closed orbits identity");

            ASSERT_EQF(kepler_orbit_period(&elements), kepler_orbit_period(&e2),
                "Closed orbits period identity");

            double P = kepler_orbit_period(&elements);
            double dt = elements.periapsis_time - e2.periapsis_time;
            double d = dt / P;

            ASSERT(ZEROF(fmod(d+0.5, 1.0)-0.5),
                "Closed orbits periapsis time identity");
        } else {
            ASSERT(!kepler_orbit_closed(&e2),
                "Open orbits identity");

            double M1 = kepler_orbit_mean_anomaly_at_time(&e2, epoch);

            if(pass == 0) {
                double ff = kepler_anomaly_mean_to_true(e2.eccentricity, M1);
                ASSERT_EQF(ff, f1,
                    "Open orbits periapsis time identity");
            } else {
                double EE = kepler_anomaly_mean_to_eccentric(e2.eccentricity, M1);
                ASSERT_EQF(EE, E1,
                    "Open orbits periapsis time identity");
            }

        }
    }

    double pos[3], vel[3], acc[3];
    kepler_elements_to_state_t(&elements, t1, pos, vel);
    double r = mag(pos);
    for(int i = 0; i < 3; ++i)
        acc[i] = - mu/(r*r*r) * pos[i];

    double dt = 1.0e-3 * M_PI/180.0 / n;
    double posplus[3], velplus[3], posminus[3], velminus[3];
    kepler_elements_to_state_t(&elements, t1-dt, posminus, velminus);
    kepler_elements_to_state_t(&elements, t1+dt, posplus, velplus);

    double dx[3], dv[3], delta_x[3], delta_v[3];
    for(int i = 0; i < 3; ++i) {
        dx[i] = 2.0 * vel[i] * dt;
        dv[i] = 2.0 * acc[i] * dt;

        delta_x[i] = posplus[i] - posminus[i];
        delta_v[i] = velplus[i] - velminus[i];
    }

    ASSERT(eqv3(dx, delta_x),
        "v = dx/dt");
    ASSERT(eqv3(dv, delta_v),
        "a = dv/dt");
}
