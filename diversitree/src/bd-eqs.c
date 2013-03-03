/*
  Constant rate Birth-death functions
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "util-splines.h"

/* For CVODES */
#include <nvector/nvector_serial.h>
#include <user_data.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_bd(double *pars, double *y, double *ydot) {
  double E = y[0], D = y[1];
  double lambda = pars[0], mu = pars[1];

  ydot[0] = mu - (mu + lambda)*E +   lambda*E*E;
  ydot[1] =    - (mu + lambda)*D + 2*lambda*D*E;
}

/* Plain BD */
/* deSolve / LSODA */
static double parms_bd[2];

void initmod_bd(void (* odeparms)(int *, double *)) {
  int N = 2;
  odeparms(&N, parms_bd);
}

void derivs_bd(int *neq, double *t, double *y, double *ydot, 
	       double *yout, int *ip) {
  do_derivs_bd(parms_bd, y, ydot);
}

/* CVODES */
int derivs_bd_cvode(realtype t, N_Vector y, N_Vector ydot,
		    void *user_data) {
  do_derivs_bd(((UserData*) user_data)->p,
	       NV_DATA_S(y),
	       NV_DATA_S(ydot));
  return 0;
}

void initial_conditions_bd(int neq, double *vars_l, double *vars_r,
			   double *pars, double t, 
			   double *vars_out) {
  vars_out[0] = vars_l[0];
  vars_out[1] = vars_l[1] * vars_r[1] * pars[0];
}
