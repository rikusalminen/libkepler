#include <libkepler/universal.h>

#include "numtest.h"

void universal_test(double *params, int num_params, void *extra_args, struct numtest_ctx *test_ctx) {
    (void)extra_args;
    ASSERT(num_params == 1, "");

    (void)params;
}
