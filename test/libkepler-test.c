#include "numtest.h"

extern numtest_callback
    anomaly_test,
    orbit_from_state_test,
    orbit_from_elements_test,
    true_anomaly_test;

const struct numtest_case numtest_cases[] = {
    { "anomaly", anomaly_test, 2, 0 },
    { "orbit_from_state", orbit_from_state_test, 7, 0 },
    { "orbit_from_elements", orbit_from_elements_test, 4, 0 },
    { "true_anomaly", true_anomaly_test, 4, 0 },
    { 0, 0, 0, 0 }
};

int main(int argc, char *argv[]) { return numtest_main(argc, argv); }
