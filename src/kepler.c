#include <stdbool.h>
#include <math.h>
#include <float.h>

#include <libkepler/kepler.h>
#include <libkepler/anomaly.h>

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
static double cube(double x) { return x*x*x; }

static void transpose3x3(double *m)
{
    for(int i = 0; i < 2; ++i) {
        for(int j = i; j < 3; ++j) {
            double temp = m[i * 3 + j];
            m[i * 3 + j] = m[j * 3 + i];
            m[j * 3 + i] = temp;
        }
    }
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

bool kepler_orbit_parabolic(const struct kepler_elements *elements) {
    return zero(elements->eccentricity - 1.0);
}

bool kepler_orbit_hyperbolic(const struct kepler_elements *elements) {
    return elements->eccentricity > 1.0 && !kepler_orbit_parabolic(elements);
}

bool kepler_orbit_closed(const struct kepler_elements *elements) {
    return elements->eccentricity < 1.0 && !kepler_orbit_parabolic(elements);
}

bool kepler_orbit_circular(const struct kepler_elements *elements) {
    return zero(elements->eccentricity);
}

double kepler_orbit_semi_latus_rectum(const struct kepler_elements *elements) {
    return elements->semi_latus_rectum;
}

double kepler_orbit_eccentricity(const struct kepler_elements *elements) {
    return elements->eccentricity;
}

double kepler_orbit_mean_motion(const struct kepler_elements *elements) {
    return elements->mean_motion;
}

double kepler_orbit_inclination(const struct kepler_elements *elements) {
    return elements->inclination;
}

double kepler_orbit_longitude_of_ascending_node(const struct kepler_elements *elements) {
    return elements->longitude_of_ascending_node;
}

double kepler_orbit_argument_of_periapsis(const struct kepler_elements *elements) {
    return elements->argument_of_periapsis;
}

double kepler_orbit_periapsis_time(const struct kepler_elements *elements) {
    return elements->periapsis_time;
}

double kepler_orbit_semi_major_axis(const struct kepler_elements *elements) {
    if(kepler_orbit_parabolic(elements))
        return INFINITY;
    return elements->semi_latus_rectum / (1.0 - square(elements->eccentricity));
}

double kepler_orbit_semi_minor_axis(const struct kepler_elements *elements) {
    if(kepler_orbit_parabolic(elements))
        return INFINITY;
    if(elements->eccentricity > 1.0)
        return -elements->semi_latus_rectum / sqrt(square(elements->eccentricity) - 1.0);
    return elements->semi_latus_rectum / sqrt(1.0 - square(elements->eccentricity));
}

double kepler_orbit_gravity_parameter(const struct kepler_elements *elements) {
    if(kepler_orbit_parabolic(elements))
        return square(elements->mean_motion) * cube(elements->semi_latus_rectum);
    return square(elements->mean_motion) *
        cube(fabs(kepler_orbit_semi_major_axis(elements)));
}

double kepler_orbit_specific_orbital_energy(const struct kepler_elements *elements) {
    if(kepler_orbit_parabolic(elements))
        return 0.0;
    return -kepler_orbit_gravity_parameter(elements) /
        (2.0 * kepler_orbit_semi_major_axis(elements));
}

double kepler_orbit_specific_angular_momentum(const struct kepler_elements *elements) {
    return sqrt(elements->semi_latus_rectum * kepler_orbit_gravity_parameter(elements));
}

double kepler_orbit_apoapsis(const struct kepler_elements *elements) {
    if(elements->eccentricity >= 1.0)
        return INFINITY;
    return elements->semi_latus_rectum / (1 - elements->eccentricity);
}

double kepler_orbit_periapsis(const struct kepler_elements *elements) {
    return elements->semi_latus_rectum / (1 + elements->eccentricity);
}

double kepler_orbit_apoapsis_vel(const struct kepler_elements *elements) {
    if(elements->eccentricity >= 1.0)
        return INFINITY;

    double mu = kepler_orbit_gravity_parameter(elements);
    return sqrt((mu / elements->semi_latus_rectum) * square(1.0 - elements->eccentricity));
}

double kepler_orbit_periapsis_vel(const struct kepler_elements *elements) {
    double mu = kepler_orbit_gravity_parameter(elements);
    return sqrt((mu / elements->semi_latus_rectum) * square(1.0 + elements->eccentricity));
}

double kepler_orbit_period(const struct kepler_elements *elements) {
    if(!kepler_orbit_closed(elements))
        return INFINITY;
    return 2.0 * M_PI / elements->mean_motion;
}

double kepler_orbit_mean_anomaly_at_time(const struct kepler_elements *elements, double t) {
    double M = elements->mean_motion * (t - elements->periapsis_time);
    if(kepler_orbit_closed(elements) && fabs(M) >= M_PI) {
        double x = (M+M_PI)/(2.0*M_PI);
        M = -M_PI + 2.0*M_PI * (x - floor(x));
    }
    return M;
}

void kepler_orientation_normal(double i, double an, double arg, double *dir) {
    (void)arg;
    dir[0] = sin(an) * sin(i);
    dir[1] = -cos(an) * sin(i);
    dir[2] = cos(i);
}

void kepler_orientation_tangent(double i, double an, double arg, double *dir) {
    dir[0] = (cos(arg) * cos(an)) - (sin(arg) * sin(an) * cos(i));
    dir[1] = (sin(arg) * cos(an) * cos(i)) + (cos(arg) * sin(an));
    dir[2] = sin(arg) * sin(i);
}

void kepler_orientation_bitangent(double i, double an, double arg, double *dir) {
    dir[0] = -(cos(arg) * sin(an) * cos(i)) - (sin(arg) * cos(an));
    dir[1] = (cos(arg) * cos(an) * cos(i)) - (sin(arg) * sin(an));
    dir[2] = cos(arg) * sin(i);
}

void kepler_orientation_matrix(double i, double an, double arg, double *mat) {
    kepler_orientation_tangent(i, an, arg, mat+0);
    kepler_orientation_bitangent(i, an, arg, mat+3);
    kepler_orientation_normal(i, an, arg, mat+6);
    transpose3x3(mat);
}

void kepler_orbit_normal(const struct kepler_elements *elements, double *dir) {
    kepler_orientation_normal(
        elements->inclination,
        elements->longitude_of_ascending_node,
        elements->argument_of_periapsis,
        dir);
}

void kepler_orbit_tangent(const struct kepler_elements *elements, double *dir) {
    kepler_orientation_tangent(
        elements->inclination,
        elements->longitude_of_ascending_node,
        elements->argument_of_periapsis,
        dir);
}

void kepler_orbit_bitangent(const struct kepler_elements *elements, double *dir) {
    kepler_orientation_bitangent(
        elements->inclination,
        elements->longitude_of_ascending_node,
        elements->argument_of_periapsis,
        dir);
}

void kepler_orbit_matrix(const struct kepler_elements *elements, double *mat) {
    kepler_orientation_matrix(
        elements->inclination,
        elements->longitude_of_ascending_node,
        elements->argument_of_periapsis,
        mat);
}

void kepler_elements_from_state(
    double mu,
    const double *pos,
    const double *vel,
    double epoch,
    struct kepler_elements *elements) {
    double r = mag(pos);
    double v2 = dot(vel, vel);

    // specific angular momentum
    double h[3];
    cross(pos, vel, h);
    // TODO: check for radial trajectory

    // eccentricity vector, direction: to periapsis, magnitude: eccentricity
    // e = 1/mu * (v^2 - mu/r) * r - dot(r, v) * v;
    double ecc[3];
    for(int i = 0; i < 3; ++i)
        ecc[i] = (1.0 / mu) * (pos[i]*(v2 - mu/r) - vel[i]*dot(pos, vel));

    double e = mag(ecc);
    bool circular = zero(e);
    bool parabolic = zero(e - 1.0);

    // line of nodes, pointing to ascending node, equatorial -> zero
    double nodes[3] = { -h[1], h[0], 0.0 };
    double N = mag(nodes);
    bool equatorial = zero(dot(nodes, nodes));

    // semi-latus rectum
    double p = dot(h, h) / mu;

    // inclination
    double i = acos(clamp(-1.0, 1.0, h[2] / mag(h)));

    // longitude of ascending node
    double an = equatorial ?  0.0 : atan2(nodes[1], nodes[0]);

    // argument of periapsis
    double arg = 0.0 / 0.0; // NaN
    if(circular)
        // circular, zero
        arg = 0.0;
    else if(equatorial)
        // equatorial, measure from X-axis, negative for retrograde
        arg = sign(h[2]) * atan2(ecc[1], ecc[0]);
    else
        // angle between eccentricity vector and line of nodes (ascending node)
        arg = sign(ecc[2]) * acos(clamp(-1.0, 1.0, dot(nodes, ecc) / (N * e)));

    // true anomaly
    double f = 0.0 / 0.0; // NaN

    if(circular && equatorial)
        // circular, equatorial -> measure from X-axis, negative for retrograde
        f = sign(h[2]) * atan2(pos[1], pos[0]);
    else if(circular)
        // circular orbit -> measure from ascending node
        f = -sign(dot(vel, nodes)) *
            acos(clamp(-1.0, 1.0, dot(nodes, pos) / (N * r)));
    else
        // measure true anomaly from periapsis (eccentricity vector)
        f = -sign(dot(vel, ecc)) *
            acos(clamp(-1.0, 1.0, dot(ecc, pos) / (e * r)));

    // mean anomaly at epoch
    double M0 = anomaly_true_to_mean(e, f);

    // mean motion
    double a = p / (1.0 - e*e);
    double n = parabolic ?
        sqrt(mu / (p*p*p)) :
        sqrt(mu / fabs(a*a*a));

    // time at periapsis
    double periapsis_time = epoch - M0 / n;

    elements->semi_latus_rectum = p;
    elements->eccentricity = e;
    elements->mean_motion = n;
    elements->inclination = i;
    elements->longitude_of_ascending_node = an;
    elements->argument_of_periapsis = arg;
    elements->periapsis_time = periapsis_time;
}

void kepler_elements_to_state_f(
    const struct kepler_elements *elements,
    double f,
    double *pos,
    double *vel) {
    // generic conic trajectory with true anomaly
    double p = elements->semi_latus_rectum;
    double e = elements->eccentricity;
    double mu = kepler_orbit_gravity_parameter(elements);

    double r = p / (1.0 + e*cos(f));
    double x = r * cos(f);
    double y = r * sin(f);

    double vx = sqrt(mu / p) * -sin(f);
    double vy = sqrt(mu / p) * (e + cos(f));

    pos[0] = x; pos[1] = y; pos[2] = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = 0.0;
}

void kepler_elements_to_state_E(
    const struct kepler_elements *elements,
    double E,
    double *pos,
    double *vel) {
    double e = elements->eccentricity;
    double p = elements->semi_latus_rectum;
    double a = kepler_orbit_semi_major_axis(elements);
    double b = kepler_orbit_semi_minor_axis(elements);
    double n = elements->mean_motion;

    double Edot = n * anomaly_dEdM(e, E);

    double x, y, vx, vy;
    if(kepler_orbit_parabolic(elements)) {
        // parabolic trajectory
        x = p/2.0 * (1.0 - E*E);
        y = p * E;

        vx = -p * E * Edot;
        vy = p * Edot;

    } else if(kepler_orbit_closed(elements)) {
        // elliptic trajectory
        x = a * (cos(E) - e);
        y = b * sin(E);

        vx = -a * sin(E) * Edot;
        vy = b * cos(E) * Edot;
    } else {
        // hyperbolic trajectory
        x = a * (cosh(E) - e);
        y = -b * sinh(E);

        vx = a * sinh(E) * Edot;
        vy = -b * cosh(E) * Edot;
    }

    pos[0] = x; pos[1] = y; pos[2] = 0.0;
    vel[0] = vx; vel[1] = vy; vel[2] = 0.0;
}

void kepler_elements_to_state_t(
    const struct kepler_elements *elements,
    double t,
    double *pos,
    double *vel) {
    // mean anomaly
    double M = kepler_orbit_mean_anomaly_at_time(elements, t);

    // eccentric anomaly
    double E = anomaly_mean_to_eccentric(elements->eccentricity, M);

    // position and velocity in orbital plane
    double pos2d[3], vel2d[3];
    kepler_elements_to_state_E(elements, E, pos2d, vel2d);

    // position and velocity in 3d
    double matrix[9];
    kepler_orbit_matrix(elements, matrix);
    matrix_vector_product(matrix, pos2d, pos);
    matrix_vector_product(matrix, vel2d, vel);
}
