#include <R.h>

static double parms_bissenew[6];

void initmod_bissenew(void (* odeparms)(int *, double *)) {
  int N = 6;
  odeparms(&N, parms_bissenew);
}

void derivs_bissenew(int *neq, double *t, double *y, double *ydot, 
                    double *yout, int *ip) {
  double E0 = y[0], E1 = y[1];
  double D0 = y[2], D1 = y[3];

  double la0 = parms_bissenew[0], la1 = parms_bissenew[1], 
    mu0 = parms_bissenew[2], mu1 = parms_bissenew[3],
    q01 = parms_bissenew[4], q10 = parms_bissenew[5];

  ydot[0] = -(mu0 + q01 + la0) * E0 + la0 * E0 * E0 + mu0 + q01 * E1;
  ydot[1] = -(mu1 + q10 + la1) * E1 + la1 * E1 * E1 + mu1 + q10 * E0;
  ydot[2] = -(mu0 + q01 + la0) * D0 + 2 * la0 * E0 * D0 + q01 * D1;
  ydot[3] = -(mu1 + q10 + la1) * D1 + 2 * la1 * E1 * D1 + q10 * D0;
}
