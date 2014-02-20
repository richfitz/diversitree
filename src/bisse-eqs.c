/*
 * These are the BiSSE equations, implemented in c
 */
#include <R.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_bisse(double *pars, const double *y, double *ydot) {
  double E0 = y[0], E1 = y[1], D0 = y[2], D1 = y[3];
  double la0 = pars[0], la1 = pars[1], mu0 = pars[2], mu1=pars[3],
    q01 = pars[4], q10 = pars[5];

  ydot[0] = -(mu0 + q01 + la0) * E0 + la0 * E0 * E0 + mu0 + q01 * E1;
  ydot[1] = -(mu1 + q10 + la1) * E1 + la1 * E1 * E1 + mu1 + q10 * E0;
  ydot[2] = -(mu0 + q01 + la0) * D0 + 2 * la0 * E0 * D0 + q01 * D1;
  ydot[3] = -(mu1 + q10 + la1) * D1 + 2 * la1 * E1 * D1 + q10 * D0;
}

/* For split models, compute some auxilliary variables  (just the
   extinction rates) */
void do_derivs_bisse_aux(double *pars, const double *y, double *ydot) {
  double E0 = y[0], E1 = y[1];
  const double la0 = pars[0], la1 = pars[1], 
    mu0 = pars[2], mu1=pars[3],
    q01 = pars[4], q10 = pars[5];

  ydot[0] = -(mu0 + q01 + la0) * E0 + la0 * E0 * E0 + mu0 + q01 * E1;
  ydot[1] = -(mu1 + q10 + la1) * E1 + la1 * E1 * E1 + mu1 + q10 * E0;
}

/* 
   Wrap these up for GslOde, which allows time and neq to be here too.
   
   TODO: Eventually move the contents of do_derivs_* into these
   functions, as no other backend possible now.
 */
void derivs_bisse_gslode(int neqs, double t, double *pars, 
			 const double *y, double *dydt) {
  do_derivs_bisse(pars, y, dydt);
}

void derivs_bisse_aux_gslode(int neqs, double t, double *pars, 
			     const double *y, double *dydt) {
  do_derivs_bisse_aux(pars, y, dydt);
}

/* 
   Not used yet, but initial condition calculations 
   
   TODO: Do we want 't' in here at all?
 */
void initial_conditions_bisse(int neq, double *vars_l, double *vars_r,
			      double *pars, double t, 
			      double *vars_out) {
  vars_out[0] = vars_l[0];
  vars_out[1] = vars_l[1];
  vars_out[2] = vars_l[2] * vars_r[2] * pars[0];
  vars_out[3] = vars_l[3] * vars_r[3] * pars[1];
}
