/* Multi-state BiSSE (MuSSE) compiled code */
#include <R.h>
#include <R_ext/BLAS.h>
#include "util.h"

/* This is the core function that actually evaluates the deriative */
void do_derivs_musse(int k, double *pars, const double *y, double *ydot) {
  const double *e = y, *d = y + k;
  double *dEdt = ydot, *dDdt = ydot + k;
  double *lambda = pars, *mu = pars + k, *Q = pars + 2*k;
  double tmp, ei, di;
  int i;

  for ( i = 0; i < k; i++ ) {
    ei = e[i];
    di = d[i];
    tmp = - lambda[i] - mu[i];
    dEdt[i] = mu[i] + tmp * ei +     lambda[i] * ei * ei;
    dDdt[i] =         tmp * di + 2 * lambda[i] * ei * di;
  }

  /* TODO: Replace with mult_mv() from util-matrix.c */
  /* mult_mv(k, Q, y, 1.0, ydot); */
  do_gemm2(Q, k, k, y, k, 2, ydot);
}

void do_derivs_musse_aux(int k, double *pars, const double *y, double *ydot) {
  const double *e = y;
  double *dEdt = ydot, ei;
  double *lambda = pars, *mu = pars + k, *Q = pars + 2*k;
  int i;

  for ( i = 0; i < k; i++ ) {
    ei = e[i];
    dEdt[i] = mu[i] - (lambda[i] + mu[i]) * ei + lambda[i] * ei * ei;
  }

  do_gemm2(Q, k, k, y, k, 1, ydot);
}

/* 
   Wrap these up for GslOde, which allows time and neq to be here too.
   
   TODO: Eventually move the contents of do_derivs_* into these
   functions, as no other backend possible now.
 */
void derivs_musse_gslode(int neqs, double t, double *pars, 
			 const double *y, double *dydt) {
  do_derivs_musse(neqs/2, pars, y, dydt);
}

void derivs_musse_aux_gslode(int neqs, double t, double *pars, 
			     const double *y, double *dydt) {
  do_derivs_musse_aux(neqs, pars, y, dydt);
}

void initial_conditions_musse(int neq, double *vars_l, double *vars_r,
			      double *pars, double t, 
			      double *vars_out) {
  const int k = neq/2;
  int i, j;
  /* E: */
  memcpy(vars_out, vars_l, k * sizeof(double));
  /* D: */
  for ( i = 0, j = k; i < k; i++, j++ )
    vars_out[j] = vars_l[j] * vars_r[j] * pars[i];
}
