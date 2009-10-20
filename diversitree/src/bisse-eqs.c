/*
 * These are the BiSSE equations, implemented in c
 */

#include <R.h>
static double parms[6];
#define lambda0 parms[0]
#define lambda1 parms[1]
#define mu0 parms[2]
#define mu1 parms[3]
#define q01 parms[4]
#define q10 parms[5]

void initmod(void (* odeparms)(int *, double *)) {
  int N = 6;
  odeparms(&N, parms);
}

void derivs(int *neq, double *t, double *y, double *ydot, 
	    double *yout, int *ip) {
  double E0 = y[0];
  double E1 = y[1];
  double D0 = y[2];
  double D1 = y[3];

  ydot[0] = -(mu0 + q01 + lambda0) * E0 + lambda0 * E0 * E0 + mu0 + q01 * E1;
  ydot[1] = -(mu1 + q10 + lambda1) * E1 + lambda1 * E1 * E1 + mu1 + q10 * E0;
  ydot[2] = -(mu0 + q01 + lambda0) * D0 + 2 * lambda0 * E0 * D0 + q01 * D1;
  ydot[3] = -(mu1 + q10 + lambda1) * D1 + 2 * lambda1 * E1 * D1 + q10 * D0;
}

void jac(int *neq, double *t, double *y, int *ml, int *mu, 
	 double *pd, int *nrowpd, double *yout, int *ip) {
  double E0 = y[0];
  double E1 = y[1];
  double D0 = y[2];
  double D1 = y[3];
  int n = *nrowpd;

  pd[0] = -(mu0 + q01 + lambda0) + 2 * lambda0 * E0;
  pd[1] = q01;

  pd[n] = q10;
  pd[n+1] = -(mu1 + q10 + lambda1) + 2 * lambda1 * E1;

  pd[2*n] = 2 * D0 * lambda0;
  pd[2*n+2] = -(mu0 + q01 + lambda0) + 2 * lambda0 * E0;
  pd[2*n+3] = q01;

  pd[3*n+1] = 2 * D1 * lambda1;
  pd[3*n+2] = q01;
  pd[3*n+3] = -(mu1 + q10 + lambda1) + 2 * lambda1 * E1;
}

