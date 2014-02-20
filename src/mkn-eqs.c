/* Mkn [Ode] - Markov k-state n-parameter character model: ODE version */
#include <R.h>
#include "util.h"

/* This is the core function that actually evaluates the deriative */
void do_derivs_mknode(int k, double *pars, const double *y, double *ydot) {
  double *Q = pars;

  /* TODO: Replace with mult_mv() from util-matrix.c */
  /* mult_mv(k, Q, y, 1.0, ydot); */
  do_gemm(Q, k, k, y, k, 1, ydot);
}

void derivs_mknode_gslode(int neqs, double t, double *pars, 
			  const double *y, double *dydt) {
  const int k = neqs;
  do_derivs_mknode(k, pars, y, dydt);
}

void initial_conditions_mknode(int neq, double *vars_l, double *vars_r,
			       double *pars, double t, 
			       double *vars_out) {
  const int k = neq;
  int i;
  for ( i = 0; i < k; i++ )
    vars_out[i] = vars_l[i] * vars_r[i];
}
