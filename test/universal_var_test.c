#include <math.h>

#include <libkepler/universal_var.h>

#include "numtest.h"

void universal_var_CS_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 1, "");

    const double maxz = 4.0*M_PI;
    double z = (-1.0 + 2.0 * params[0]) * maxz;
    double cf = universal_var_C(z), cs = universal_var_Cseries(z);
    double sf = universal_var_S(z), ss = universal_var_Sseries(z);

    ASSERT(isfinite(cf) && isfinite(cs) && isfinite(sf) && isfinite(ss),
        "C and S not NaN");

    ASSERT_EQF(cs, cf, "C series equals C function w/trigonometry");
    ASSERT_EQF(ss, sf, "S series equals S function w/trigonometry");

    double dCdz = universal_var_dCdz(z, cs, ss);
    double dSdz = universal_var_dSdz(z, cs, ss);

    if(ZEROF(z)) // TODO: derivatives at z = 0
        dCdz = dSdz = 0.0;

    ASSERT(isfinite(dCdz) && isfinite(dSdz),
        "dC/dz and dS/dz not NaN");

    double dz = 1.0e-10 * maxz;
    double Cplus = universal_var_Cseries(z + dz);
    double Cminus = universal_var_Cseries(z - dz);
    double Splus = universal_var_Sseries(z + dz);
    double Sminus = universal_var_Sseries(z - dz);

    ASSERT_EQF(2.0 * dCdz * dz, Cplus - Cminus, "dC/dz");
    ASSERT_EQF(2.0 * dSdz * dz, Splus - Sminus, "dS/dz");
}
