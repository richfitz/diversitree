/* Either BBM or meristic states with just an up and down rate */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include "util.h"

/* For CVODES */
#include <nvector/nvector_serial.h>
#include <user_data.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_mkn_meristic(int k, double *pars, double *y, double *ydot) {
  const double d = pars[0], u = pars[1];
  const double ud = u + d;
  double uy, dy;
  int i;

  uy = u * y[0];
  dy = d * y[1];
  
  ydot[0] = -u * y[0] + u * y[1];
  for ( i = 1; i < k - 1; i++ )
    ydot[i] = -ud * y[i] + d * y[i-1] + u * y[i+1];
  ydot[k-1] = -d * y[k-1] + d * y[k-2];
}

/* deSolve / LSODA */
static double *parms_mkn_meristic;

void initmod_mkn_meristic(void (* odeparms)(int *, double *)) {
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  parms_mkn_meristic = REAL(get_deSolve_gparms());
} 

void derivs_mkn_meristic(int *neq, double *t, double *y, double *ydot,
			 double *yout, int *ip) {
  do_derivs_mkn_meristic(*neq, parms_mkn_meristic, y, ydot);
}

/* CVODES */
int derivs_mkn_meristic_cvode(realtype t, N_Vector y, N_Vector ydot,
			      void *user_data) {
  const UserData *data = (UserData*) user_data;
  do_derivs_mkn_meristic(data->neq,
			 data->p,
			 NV_DATA_S(y),
			 NV_DATA_S(ydot));
  return 0;
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
