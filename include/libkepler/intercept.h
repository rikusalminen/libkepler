#ifndef LIBKEPLER_INTERCEPT_H
#define LIBKEPLER_INTERCEPT_H

#include <stdbool.h>

struct kepler_elements;

struct intercept {
    double time;
    double distance;
    double position[3];
    double relative_velocity[3];
};

bool intercept_orbit(
    const struct kepler_elements *orbit1,
    const struct kepler_elements *orbit2,
    double t0, double t1,
    struct intercept *intercept);

#endif
