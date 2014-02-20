/* Either BBM or meristic states with just an up and down rate */
#include <R.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_mkn_meristic(int k, double *pars, const double *y, 
			    double *ydot) {
  const double d = pars[0], u = pars[1];
  const double ud = u + d;
  int i;

  ydot[0] = -u * y[0] + u * y[1];
  for ( i = 1; i < k - 1; i++ )
    ydot[i] = -ud * y[i] + d * y[i-1] + u * y[i+1];
  ydot[k-1] = -d * y[k-1] + d * y[k-2];
}

void derivs_mkn_meristic_gslode(int neqs, double t, double *pars, 
				const double *y, double *dydt) {
  const int k = neqs;
  do_derivs_mkn_meristic(k, pars, y, dydt);
}

void initial_conditions_mkn_meristic(int neq, 
				     double *vars_l, double *vars_r,
				     double *pars, double t, 
				     double *vars_out) {
  const int k = neq;
  int i;
  for ( i = 0; i < k; i++ )
    vars_out[i] = vars_l[i] * vars_r[i];
}
