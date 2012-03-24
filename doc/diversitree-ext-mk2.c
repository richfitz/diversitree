#include <R.h>

static double parms_mk2new[2];

void initmod_mk2new(void (* odeparms)(int *, double *)) {
  int N = 2;
  odeparms(&N, parms_mk2new);
}

void derivs_mk2new(int *neq, double *t, double *y, double *ydot, 
                    double *yout, int *ip) {
  double D0 = y[0], D1 = y[1];
  double q01 = parms_mk2new[0], q10 = parms_mk2new[1];

  ydot[0] = -q01 * D0 + q01 * D1;
  ydot[1] = -q10 * D1 + q10 * D0;
}
