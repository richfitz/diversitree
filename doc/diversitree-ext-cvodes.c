#include <R.h>
/* For CVODES */
#include <nvector/nvector_serial.h>
#include <user_data.h>

/* This is the core function that actually evaluates the deriative, as
   in the previous version.  However, to save on duplication, I've
   separated this out a little. */
void do_derivs_bissenew(double *pars, double *y, double *ydot) {
  double E0 = y[0], E1 = y[1], D0 = y[2], D1 = y[3];
  double la0 = pars[0], la1 = pars[1], mu0 = pars[2], mu1=pars[3],
    q01 = pars[4], q10 = pars[5];

  ydot[0] = -(mu0 + q01 + la0) * E0 + la0 * E0 * E0 + mu0 + q01 * E1;
  ydot[1] = -(mu1 + q10 + la1) * E1 + la1 * E1 * E1 + mu1 + q10 * E0;
  ydot[2] = -(mu0 + q01 + la0) * D0 + 2 * la0 * E0 * D0 + q01 * D1;
  ydot[3] = -(mu1 + q10 + la1) * D1 + 2 * la1 * E1 * D1 + q10 * D0;
}

/* Here is the deSolve interface, as before, but using the
   do_derivs_bissenew function defined above */
static double parms_bissenew[6];
void initmod_bissenew(void (* odeparms)(int *, double *)) {
  int N = 6;
  odeparms(&N, parms_bissenew);
}
void derivs_bissenew(int *neq, double *t, double *y, double *ydot, 
		     double *yout, int *ip) {
  do_derivs_bissenew(parms_bissenew, y, ydot);
}

/* CVODES */
int derivs_bissenew_cvode(realtype t, N_Vector y, N_Vector ydot,
			  void *user_data) {
  do_derivs_bissenew(((UserData*) user_data)->p,
		     NV_DATA_S(y),
		     NV_DATA_S(ydot));
  return 0;
}

/* This is also required, to compute initial conditions.  It is almost
   identical to the R version */
void initial_conditions_bissenew(int neq, double *vars_l, double *vars_r,
				 double *pars, double t, 
				 double *vars_out) {
  vars_out[0] = vars_l[0]; /* E0, first branch only */
  vars_out[1] = vars_l[1]; /* E1, first branch only */
  vars_out[2] = vars_l[2]*vars_r[2]*pars[0]; /* D0_l*D0_r*lambda0 */
  vars_out[3] = vars_l[3]*vars_r[3]*pars[1]; /* D1_l*D1_r*lambda1 */
}
