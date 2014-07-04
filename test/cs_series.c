#include <math.h>
#include <stdio.h>

#include <libkepler/universal_var.h>

int main(int argc, char *argv[]) {
    unsigned num_steps = 0;

    if(argc < 2 || sscanf(argv[1], "%u\n", &num_steps) != 1)
        num_steps = 101;

    printf("# z:\tC:\tS:\tdC/dz:\tdS/dz:\n");
    for(unsigned i = 0; i < num_steps; ++i) {
        double t = -1.0 + 2.0 * (i / (double)(num_steps-1));
        double maxz = 16.0 * M_PI * M_PI;
        double z = t * maxz;

        double C = universal_var_Cseries(z);
        double S = universal_var_Sseries(z);
        double dCdz = universal_var_dCdz(z, C, S);
        double dSdz = universal_var_dSdz(z, C, S);

        printf("%lf\t%lf\t%lf\t%lf\t%lf\n", z, C, S, dCdz, dSdz);
    }

    return 0;
}
