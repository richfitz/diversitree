/* MknOde - Markov k-state n-parameter character model: ODE version */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include "util.h"

/* For CVODES */
#include <nvector/nvector_serial.h>
#include <user_data.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_mknode(int k, double *pars, double *y, double *ydot) {
  double *Q = pars;

  /* TODO: Replace with mult_mv() from util-matrix.c */
  /* mult_mv(k, Q, y, 1.0, ydot); */
  do_gemm(Q, k, k, y, k, 1, ydot);
}

/* Plain Mk */
/* deSolve / LSODA */
static double *parms_mknode;

void initmod_mknode(void (* odeparms)(int *, double *)) {
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  parms_mknode = REAL(get_deSolve_gparms());
} 

void derivs_mknode(int *neq, double *t, double *y, double *ydot,
		double *yout, int *ip) {
  do_derivs_mknode(*neq, parms_mknode, y, ydot);
}

/* CVODES */
int derivs_mknode_cvode(realtype t, N_Vector y, N_Vector ydot,
		     void *user_data) {
  const UserData *data = (UserData*) user_data;
  do_derivs_mknode(data->neq,
		data->p,
		NV_DATA_S(y),
		NV_DATA_S(ydot));
  return 0;
}

void initial_conditions_mknode(int neq, double *vars_l, double *vars_r,
			    double *pars, double t, 
			    double *vars_out) {
  const int k = neq;
  int i;
  for ( i = 0; i < k; i++ )
    vars_out[i] = vars_l[i] * vars_r[i];
}
