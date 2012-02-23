/*
 * 
 *
 * bisse_ness-eqs.c
 *  
 *
 * Created by Karen Magnuson-Ford on February 16, 2012.
 * Copyright 2012. All rights reserved.
 *
 * Modified by RGF 22 Feb 2012.
 *
 */

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* For CVODES */
#include <nvector/nvector_serial.h>
#include <user_data.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_bisseness(double *pars, double *y, double *ydot) {
  double E0 = y[0], E1 = y[1], D0 = y[2], D1 = y[3];
  double la0 = pars[0], la1 = pars[1], mu0 = pars[2], mu1 = pars[3],
    q01 = pars[4], q10 = pars[5],
    p0c = pars[6], p0a = pars[7], p1c = pars[8], p1a = pars[9];

  ydot[0] = -(mu0 + q01 + la0) * E0 + 
    la0 * E0 * E0 * (1 - p0c) + mu0 + q01 * E1 + 
    la0 * E0 * E1 * p0c * p0a + la0 * E1 * E1 * p0c * (1 - p0a);
  ydot[1] = -(mu1 + q10 + la1) * E1 + 
    la1 * E1 * E1 * (1 - p1c) + mu1 + q10 * E0 + 
    la1 * E0 * E1 * p1c * p1a + la1 * E0 * E0 * p1c * (1 - p1a);
  ydot[2] = -(mu0 + q01 + la0) * D0 + 
    2 * la0 * E0 * D0 * (1 - p0c) + q01 * D1 + 
    (E0 * D1 + E1 * D0) * la0 * p0c * p0a + 
    2 * la0 * E1 * D1 * p0c * (1 - p0a);
  ydot[3] = -(mu1 + q10 + la1) * D1 + 
    2 * la1 * E1 * D1 * (1 - p1c) + q10 * D0 + 
    (E1 * D0 + E0 * D1) * la1 * p1c * p1a + 
    2 * la1 * E0 * D0 * p1c * (1 - p1a);
}

/* deSolve / LSODA */
static double parms_bisseness[10];

void initmod_bisseness(void (* odeparms)(int *, double *)) {
  int N = 10;
  odeparms(&N, parms_bisseness);
}

void derivs_bisseness(int *neq, double *t, double *y, double *ydot, 
		  double *yout, int *ip) {
  do_derivs_bisseness(parms_bisseness, y, ydot);
}

/* CVODES */
int derivs_bisseness_cvode(realtype t, N_Vector y, N_Vector ydot,
			   void *user_data) {
  do_derivs_bisseness(((UserData*) user_data)->p,
		      NV_DATA_S(y),
		      NV_DATA_S(ydot));
  return 0;
}

/* This actually needs completely updating */
void initial_conditions_bisseness(int neq, double *vars_l, double *vars_r,
				  double *pars, double t, 
				  double *vars_out) {
  error("this does not work yet");
  vars_out[0] = vars_l[0];
  vars_out[1] = vars_l[1];
  vars_out[2] = vars_l[2] * vars_r[2] * pars[0];
  vars_out[3] = vars_l[3] * vars_r[3] * pars[1];
}
