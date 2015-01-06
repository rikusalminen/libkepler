#include <libkepler/eccentric_anomaly.h>

#include <math.h>
#include <float.h>

#include "numtest.h"

static int zero(double x) {
    return x*x < DBL_EPSILON;
}


void eccentric_anomaly_test(
    double *params, int num_params,
    void *extra_args,
    struct numtest_ctx *test_ctx) {
    (void)extra_args;

    ASSERT(num_params == 4, "");

    double mu = 1.0 + params[0] * 1.0e10;
    double p = 1.0 + params[1] * 1.0e10;
    double e = params[2] * 4.0;

    double E = (-1.0 + params[3] * 2.0) * M_PI;

    double a = p / (1.0 - e*e);
    double b = e < 1.0 ? p / sqrt(1.0 - e*e) : -p / sqrt(e*e - 1.0);
    double c = e < 1.0 ? a*e : -a*e;
    double q = p / (1.0 + e);

    double n = zero(e - 1.0) ? sqrt(mu / (p*p*p)) : sqrt(mu / fabs(a*a*a));
    double dt = (2.0 * M_PI / n) * (1.0 / 36000.0); // 1/100 degree

    double r = eccentric_radius(p, e, E);
    ASSERT(isfinite(r), "Radius not NaN");
    ASSERT(r > 0, "Radius is positive");

    ASSERT_LTF(q, r, "Radius is larger than periapsis");
    if(e < 1.0)
        ASSERT_LTF(r, p / (1.0 - e), "Radius is smaller than apoapsis");

    if(!zero(e)) { // not a circular orbit
        double EE = eccentric_anomaly_from_radius(p, e, r);
        ASSERT(isfinite(EE), "Eccentric anomaly not NaN");
        ASSERT_RANGEF(EE, 0, M_PI, "Eccentric anomaly range");

        if(zero(E)) // Accuracy is bad near periapsis
            ASSERT(zero(EE*EE), "Eccentric anomaly radius identity (zero)");
        else
            ASSERT_EQF(fabs(E), EE, "Eccentric anomaly radius identity");
    }

    double dEdt = eccentric_dEdt(mu, p, e, E);
    double dE = dEdt * dt;
    ASSERT(isfinite(dEdt), "dE/dt not NaN");
    ASSERT(dEdt > 0.0, "dE/dt is positive");

    double t = eccentric_time(mu, p, e, E);
    ASSERT(isfinite(t), "Eccentric anomaly time of flight not NaN");

    double tplus = eccentric_time(mu, p, e, E+dE);
    double tminus = eccentric_time(mu, p, e, E-dE);
    ASSERT_EQF(tplus - tminus, 2.0*dt, "dE/dt");

    double v = eccentric_velocity(mu, p, e, E);
    double rdot = eccentric_velocity_radial(mu, p, e, E);
    double rfdot = eccentric_velocity_horizontal(mu, p, e, E);
    ASSERT(isfinite(v) && isfinite(rdot) && isfinite(rfdot),
        "Velocity not NaN");
    ASSERT_EQF(rdot*rdot + rfdot*rfdot, v*v,
        "Velocity magitude");
    ASSERT(rfdot > 0, "Horizontal velocity is positive");

    double h = sqrt(mu * p);
    ASSERT_EQF(h, r * rfdot,
        "Specific relative angular momentum");

    double energy = v*v/2.0 - mu/r;
    if(zero(e-1.0)) { // Accuracy is bad for parabolic trajectory
        ASSERT(zero(energy*energy),
            "Specific orbital energy (vis-viva) is zero (parabolic)");
    } else {
        double visviva = zero(e - 1.0) ? 0.0 : -mu/(2.0*a);
        ASSERT_EQF(visviva, energy,
            "Specific orbital energy (vis-viva)");
    }

    double rplus = eccentric_radius(p, e, E+dE);
    double rminus = eccentric_radius(p, e, E-dE);
    ASSERT_EQF(rplus - rminus, 2.0*rdot*dt,
        "rdot = dr/dt");

    double tan_phi = eccentric_tan_phi(e, E);
    double phi = eccentric_flight_path_angle(e, E);

    ASSERT(isfinite(tan_phi) && isfinite(phi),
        "Flight path angle not NaN");
    ASSERT_EQF(rdot / rfdot, tan_phi,
        "Flight path angle");
    ASSERT_EQF(tan(phi), tan_phi,
        "Flight path angle tangent");

    double x = eccentric_x(p, e, E);
    double y = eccentric_y(p, e, E);
    double xdot = eccentric_xdot(mu, p, e, E);
    double ydot = eccentric_ydot(mu, p, e, E);
    ASSERT(isfinite(x) && isfinite(y) &&
        isfinite(xdot) && isfinite(ydot),
        "Position and velocity not NaN");

    ASSERT_EQF(x*x + y*y, r*r,
        "Position magnitude");
    ASSERT_EQF(xdot*xdot + ydot*ydot, v*v,
        "Velocity magnitude (xy)");

    ASSERT_EQF(h, x * ydot - y * xdot,
        "Specific relative angular momentum (xy)");

    if(!zero(e))
        ASSERT_EQF(p/e, x + r/e,
            "Focus-Directrix property");

    if(!zero(e-1.0)) { // TODO: focal-radii property for parabolic trajectory
        double x2 = e < 1.0 ? x+2.0*c : x-2.0*c;
        double r2 = sqrt(x2*x2 + y*y);
        double sum = e < 1.0 ? r2 + r : r2 - r;

        ASSERT_EQF(e < 1.0 ? 2.0*a : -2.0*a, sum,
            "Focal-radii property");
    }

    if(zero(e-1.0)) {
        ASSERT_EQF(2.0*p*(x - q), -(y*y),
            "Trajectory is a parabola");
    } else if(e > 1.0) {
        ASSERT_EQF((x-c)*(x-c)/(a*a), 1.0 + y*y/(b*b),
            "Trajectory is a hyperbola");
    } else {
        ASSERT_EQF((x+c)*(x+c)/(a*a), 1.0 - y*y/(b*b),
            "Trajectory is an ellipse");
    }

    double dotrv = x*xdot + y*ydot;
    if(zero(e-1.0)) {
        ASSERT_EQF(dotrv / sqrt(mu*p), E,
            "Parabolic anomaly r dot v");
        ASSERT_EQF(r/q - 1.0, E*E,
            "Parabolic anomaly squared");
    } else if(e > 1.0) {
        ASSERT_EQF(dotrv / sqrt(-mu*a), e*sinh(E),
            "Hyperbolic anomaly sine");
        ASSERT_EQF(1.0 - r/a, e*cosh(E),
            "Hyperbolic anomaly cosine");
    } else {
        ASSERT_EQF(dotrv / sqrt(mu*a), e*sin(E),
            "Eccentric anomaly sine");
        ASSERT_EQF(1.0 - r/a, e*cos(E),
            "Eccentric anomaly cosine");
    }

    double xplus = eccentric_x(p, e, E+dE);
    double yplus = eccentric_y(p, e, E+dE);
    double xminus = eccentric_x(p, e, E-dE);
    double yminus = eccentric_y(p, e, E-dE);
    double dx = xplus - xminus, dy = yplus - yminus;
    double dx2 = 2.0 * xdot * dt, dy2 = 2.0 * ydot * dt;
    double xx = dx - dx2, yy = dy - dy2;

    ASSERT(ZEROF((xx*xx + yy*yy) / (dx2*dx2 + dy2*dy2)),
        "v = dr/dt");

    double acc = mu / (r*r);
    double dvx2 = -(acc * x / r) * 2.0*dt, dvy2 = -(acc * y / r) * 2.0*dt;
    double xdotplus = eccentric_xdot(mu, p, e, E+dE);
    double ydotplus = eccentric_ydot(mu, p, e, E+dE);
    double xdotminus = eccentric_xdot(mu, p, e, E-dE);
    double ydotminus = eccentric_ydot(mu, p, e, E-dE);
    double dvx = xdotplus - xdotminus, dvy = ydotplus - ydotminus;
    double vxx = dvx - dvx2, vyy = dvy - dvy2;

    ASSERT(ZEROF((vxx*vxx + vyy*vyy) / (dvx2*dvx2 + dvy2*dvy2)),
        "a = dv/dt");
}
