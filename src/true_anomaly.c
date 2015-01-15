#include <libkepler/true_anomaly.h>

#include <math.h>
#include <float.h>

static double clamp(double min, double max, double x) {
    return (x < min ? min : (x > max ? max : x));
}

static double square(double x) { return x*x; }
static double cube(double x) { return x*x*x; }


double true_radius(double p, double e, double f) {
    return p / (1.0 + e * cos(f));
}

double true_anomaly_from_radius(double p, double e, double r) {
    return acos(clamp(-1.0, 1.0, (p / r - 1.0) / e));
}

double true_dfdt(double mu, double p, double e, double f) {
    return sqrt(mu / cube(p)) * square(1 + e * cos(f));
}

double true_velocity(double mu, double p, double e, double f) {
    return sqrt(mu / p * (e*e + 2.0*e*cos(f) + 1));
}

double true_velocity_radial(double mu, double p, double e, double f) {
    return sqrt(mu/p) * e * sin(f);
}

double true_velocity_horizontal(double mu, double p, double e, double f) {
    return sqrt(mu/p) * (1 + e * cos(f));
}


double true_tan_phi(double e, double f) {
    return e * sin(f) / (1 + e * cos(f));
}

double true_flight_path_angle(double e, double f) {
    return atan(true_tan_phi(e, f));
}

double true_x(double p, double e, double f) {
    return cos(f) * true_radius(p, e, f);
}

double true_y(double p, double e, double f) {
    return sin(f) * true_radius(p, e, f);
}

double true_xdot(double mu, double p, double e, double f) {
    (void)e;
    return -sqrt(mu / p) * sin(f);
}

double true_ydot(double mu, double p, double e, double f) {
    return sqrt(mu / p) * (e + cos(f));
}
