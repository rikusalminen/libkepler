#include "numtest.h"

extern numtest_callback
    conic_test,
    anomaly_test,
    orientation_test,
    orbit_from_state_test,
    orbit_from_elements_test,
    true_anomaly_test,
    eccentric_anomaly_test,
    orbit_from_state_test_new,
    orbit_from_elements_test_new,
    dummy_test;

const struct numtest_case numtest_cases[] = {
    { "conic", conic_test, 3, 0 },
    { "anomaly", anomaly_test, 2, 0 },
    { "orientation", orientation_test, 3, 0 },
    { "orbit_from_state", orbit_from_state_test, 7, 0 },
    { "orbit_from_elements", orbit_from_elements_test, 4, 0 },
    { "true_anomaly", true_anomaly_test, 4, 0 },
    { "eccentric_anomaly", eccentric_anomaly_test, 4, 0 },
    { "orbit_from_state_new", orbit_from_state_test_new, 7, 0 },
    { "orbit_from_elements_new", orbit_from_elements_test_new, 6, 0 },
    { 0, 0, 0, 0 }
};

int main(int argc, char *argv[]) { return numtest_main(argc, argv); }
