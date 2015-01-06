#include <libkepler/eccentric_anomaly.h>

#include <math.h>
#include <float.h>

static int zero(double x) {
    return x*x < DBL_EPSILON;
}

static double clamp(double min, double max, double x) {
    return (x < min ? min : (x > max ? max : x));
}


double eccentric_radius(double p, double e, double E) {
    double a = p / (1.0 - e*e);
    double q = p / (1.0 + e);

    if(zero(e-1.0)) {
        // parabolic
        return q * (E*E + 1.0);
    } else if(e > 1.0) {
        // hyperbolic
        return a * (1.0 - e*cosh(E));
    } else {
        // elliptic
        return a * (1.0 - e*cos(E));
    }
}

double eccentric_anomaly_from_radius(double p, double e, double r) {
    double a = p / (1.0 - e*e);
    double q = p / (1.0 + e);

    if(zero(e-1.0)) {
        // parabolic
        return sqrt(r/q - 1.0);
    } else if(e > 1.0) {
        // hyperbolic
        return acosh(fmax(1.0, (1.0 - r/a) / e));
    } else {
        // elliptic
        return acos(clamp(-1.0, 1.0, (1.0 - r/a) / e));
    }
}

double eccentric_dEdt(double mu, double p, double e, double E) {
    double a = p / (1.0 - e*e);

    if(zero(e-1.0)) {
        // parabolic
        return sqrt(mu / p) / eccentric_radius(p, e, E);
    } else if(e > 1.0) {
        // hyperbolic
        return sqrt(mu / -a) / eccentric_radius(p, e, E);
    } else {
        // elliptic
        return sqrt(mu / a) / eccentric_radius(p, e, E);
    }
}

double eccentric_time(double mu, double p, double e, double E) {
    double a = p / (1.0 - e*e);

    if(zero(e-1.0)) {
        // parabolic
        return sqrt(p*p*p / mu) * (E*E*E / 6.0 + E / 2.0);
    } else if(e > 1.0) {
        // hyperbolic
        return sqrt(-a*a*a / mu) * (e*sinh(E) - E);
    } else {
        // elliptic
        return sqrt(a*a*a / mu) * (E - e*sin(E));
    }
}

double eccentric_velocity(double mu, double p, double e, double E) {
    double a = p / (1.0 - e*e);

    if(zero(e-1.0)) {
        // parabolic
        return sqrt((mu/p) * 4.0 / (E*E + 1.0));
    } else if(e > 1.0) {
        // hyperbolic
        return sqrt((mu/-a) * (e*cosh(E) + 1.0) / (e*cosh(E) - 1.0));
    } else {
        // elliptic
        return sqrt((mu/a) * (1.0 + e*cos(E)) / (1.0 - e*cos(E)));
    }
}

double eccentric_velocity_radial(double mu, double p, double e, double E) {
    double a = p / (1.0 - e*e);

    if(zero(e-1.0)) {
        // parabolic
        return sqrt(mu/p) * 2.0*E / (E*E + 1.0);
    } else if(e > 1.0) {
        // hyperbolic
        return sqrt(mu/-a) * e*sinh(E) / (e*cosh(E) - 1.0);
    } else {
        // elliptic
        return sqrt(mu/a) * e*sin(E) / (1.0 - e*cos(E));
    }
}

double eccentric_velocity_horizontal(double mu, double p, double e, double E) {
    double a = p / (1.0 - e*e);

    if(zero(e-1.0)) {
        // parabolic
        return sqrt(mu/p) * 2.0 / (E*E + 1.0);
    } else if(e > 1.0) {
        // hyperbolic
        return sqrt((mu/-a) * (e*e - 1.0)) / (e*cosh(E) - 1.0);
    } else {
        // elliptic
        return sqrt((mu/a) * (1.0 - e*e)) / (1.0 - e*cos(E));
    }
}

double eccentric_tan_phi(double e, double E) {
    if(zero(e-1.0)) {
        // parabolic
        return E;
    } else if(e > 1.0) {
        // hyperbolic
        return e*sinh(E) / sqrt(e*e - 1.0);
    } else {
        // elliptic
        return e*sin(E) / sqrt(1.0 - e*e);
    }
}

double eccentric_flight_path_angle(double e, double E) {
    return atan(eccentric_tan_phi(e, E));
}

double eccentric_x(double p, double e, double E) {
    double a = p / (1.0 - e*e);
    double q = p / (1.0 + e);

    if(zero(e-1.0)) {
        // parabolic
        return q * (1.0 - E*E);
    } else if(e > 1.0) {
        // hyperbolic
        return a * (cosh(E) - e);
    } else {
        // elliptic
        return a * (cos(E) - e);
    }
}

double eccentric_y(double p, double e, double E) {
    double b = e < 1.0 ? p/sqrt(1.0 - e*e) : p/sqrt(e*e - 1.0);

    if(zero(e-1.0)) {
        // parabolic
        return p * E;
    } else if(e > 1.0) {
        // hyperbolic
        return b * sinh(E);
    } else {
        // elliptic
        return b * sin(E);
    }
}

double eccentric_xdot(double mu, double p, double e, double E) {
    double a = p / (1.0 - e*e);

    if(zero(e-1.0)) {
        // parabolic
        return sqrt(mu / p) * -2.0*E / (E*E + 1.0);
    } else if(e > 1.0) {
        // hyperbolic
        return sqrt(mu / -(a*a*a)) * a*sinh(E) / (e*cosh(E) - 1.0);
    } else {
        // elliptic
        return sqrt(mu / (a*a*a)) * -a*sin(E) / (1.0 - e*cos(E));
    }
}

double eccentric_ydot(double mu, double p, double e, double E) {
    double a = p / (1.0 - e*e);
    double b = e < 1.0 ? p/sqrt(1.0 - e*e) : p/sqrt(e*e - 1.0);

    if(zero(e-1.0)) {
        // parabolic
        return sqrt(mu / p) * 2.0 / (E*E + 1.0);
    } else if(e > 1.0) {
        // hyperbolic
        return sqrt(mu / -(a*a*a)) * b*cosh(E) / (e*cosh(E) - 1.0);
    } else {
        // elliptic
        return sqrt(mu / (a*a*a)) * b*cos(E) / (1.0 - e*cos(E));
    }
}
