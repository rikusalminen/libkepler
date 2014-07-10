#ifndef LIBKEPLER_KEPLER_H
#define LIBKEPLER_KEPLER_H

#include <stdbool.h>

struct kepler_elements {
    double semi_latus_rectum;
    double eccentricity;
    double mean_motion;
    double inclination;
    double longitude_of_ascending_node;
    double argument_of_periapsis;
    double periapsis_time;
};

double kepler_anomaly_mean_to_eccentric(double e, double M);
double kepler_anomaly_eccentric_to_mean(double e, double E);
double kepler_anomaly_eccentric_to_true(double e, double E);
double kepler_anomaly_true_to_eccentric(double e, double f);
double kepler_anomaly_true_to_mean(double e, double f);
double kepler_anomaly_mean_to_true(double e, double M);
double kepler_anomaly_dEdM(double e, double E);

bool kepler_orbit_parabolic(const struct kepler_elements *elements);
bool kepler_orbit_hyperbolic(const struct kepler_elements *elements);
bool kepler_orbit_closed(const struct kepler_elements *elements);
bool kepler_orbit_circular(const struct kepler_elements *elements);

double kepler_orbit_semi_latus_rectum(const struct kepler_elements *elements);
double kepler_orbit_eccentricity(const struct kepler_elements *elements);
double kepler_orbit_mean_motion(const struct kepler_elements *elements);
double kepler_orbit_inclination(const struct kepler_elements *elements);
double kepler_orbit_longitude_of_ascending_node(const struct kepler_elements *elements);
double kepler_orbit_argument_of_periapsis(const struct kepler_elements *elements);
double kepler_orbit_periapsis_time(const struct kepler_elements *elements);
double kepler_orbit_semi_major_axis(const struct kepler_elements *elements);
double kepler_orbit_semi_minor_axis(const struct kepler_elements *elements);
double kepler_orbit_gravity_parameter(const struct kepler_elements *elements);
double kepler_orbit_specific_orbital_energy(const struct kepler_elements *elements);
double kepler_orbit_specific_angular_momentum(const struct kepler_elements *elements);
double kepler_orbit_apoapsis(const struct kepler_elements *elements);
double kepler_orbit_periapsis(const struct kepler_elements *elements);
double kepler_orbit_apoapsis_vel(const struct kepler_elements *elements);
double kepler_orbit_periapsis_vel(const struct kepler_elements *elements);
double kepler_orbit_period(const struct kepler_elements *elements);
double kepler_orbit_mean_anomaly_at_time(const struct kepler_elements *elements, double t);
void kepler_orientation_normal(double i, double an, double arg, double *dir);
void kepler_orientation_tangent(double i, double an, double arg, double *dir);
void kepler_orientation_bitangent(double i, double an, double arg, double *dir);
void kepler_orientation_matrix(double i, double an, double arg, double *mat);
void kepler_orbit_normal(const struct kepler_elements *elements, double *dir);
void kepler_orbit_tangent(const struct kepler_elements *elements, double *dir);
void kepler_orbit_bitangent(const struct kepler_elements *elements, double *dir);
void kepler_orbit_matrix(const struct kepler_elements *elements, double *mat);

void kepler_elements_from_state(
    double mu,
    const double *pos,
    const double *vel,
    double epoch,
    struct kepler_elements *elements);
void kepler_elements_to_state_f(
    const struct kepler_elements *elements,
    double f,
    double *pos,
    double *vel);
void kepler_elements_to_state_E(
    const struct kepler_elements *elements,
    double E,
    double *pos,
    double *vel);
void kepler_elements_to_state_t(
    const struct kepler_elements *elements,
    double t,
    double *pos,
    double *vel);

#endif
