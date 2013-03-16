/*
  Constant rate Birth-death functions
 */
#include <R.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_bd(double *pars, const double *y, double *ydot) {
  double E = y[0], D = y[1];
  double lambda = pars[0], mu = pars[1];

  ydot[0] = mu - (mu + lambda)*E +   lambda*E*E;
  ydot[1] =    - (mu + lambda)*D + 2*lambda*D*E;
}

void derivs_bd_gslode(int neqs, double t, double *pars, 
		      const double *y, double *dydt) {
  do_derivs_bd(pars, y, dydt);
}

void initial_conditions_bd(int neq, double *vars_l, double *vars_r,
			   double *pars, double t, 
			   double *vars_out) {
  vars_out[0] = vars_l[0];
  vars_out[1] = vars_l[1] * vars_r[1] * pars[0];
}
