/* Multi-state BiSSE compiled code */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include "util.h"

void do_derivs_musse(int k, double *pars, double *y, double *ydot);

static double *parms_musse;
void initmod_musse(void (* odeparms)(int *, double *)) {
  /* TODO: I should check here about the lengths of parameters, but I
     won't bother; it is not clear how best to do this, anyway.  Most
     of the checking should be done in the R end of things.  Because
     this is going to get an R object, it should be much easier to the
     checking there. */
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  parms_musse = REAL(get_deSolve_gparms());
} 

void derivs_musse(int *neq, double *t, double *y, double *ydot,
		   double *yout, int *ip) {

  do_derivs_musse(*neq / 2, parms_musse, y, ydot);
  /*
  const int k = *neq / 2;
  double *e = y, *d = y + k;
  double *dEdt = ydot, *dDdt = ydot + k;
  double *lambda = parms_musse, *mu = parms_musse + k, 
    *Q = parms_musse + 2*k;
  int i;
  double tmp, ei, di;
  for ( i = 0; i < k; i++ ) {
    ei = e[i];
    di = d[i];
    tmp = - lambda[i] - mu[i];
    dEdt[i] = mu[i] + tmp * ei +     lambda[i] * ei * ei;
    dDdt[i] =         tmp * di + 2 * lambda[i] * ei * di;
  }

  do_gemm2(Q, k, k, y, k, 2, ydot);
  */
}

/* This gives a skeleton for the final version */
void do_derivs_musse(int k, double *pars, double *y, double *ydot) {
  double *e = y, *d = y + k;
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

  do_gemm2(Q, k, k, y, k, 2, ydot);
}

/* This might go in its own file; this is the extent of the
   time-dependent MuSSE code - it's not much to look at. 

   This is based on the section in R-exts on "evaluating R expressions
   from C".
 */
static SEXP func_musse;
static SEXP rho_musse;

void initmod_musse_t(void (* odeparms)(int *, double *)) {
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  SEXP obj = get_deSolve_gparms();

  func_musse  = VECTOR_ELT(obj, 0);
  rho_musse   = VECTOR_ELT(obj, 1);

  if ( !isFunction(func_musse) )
    error("First element must be a function");
  if ( !isEnvironment(rho_musse) )
    error("Second element must be an environment");
}

void derivs_musse_t(int *neq, double *t, double *y, double *ydot, 
		    double *yout, int *ip) {
  SEXP R_fcall;
  double *pars;

  PROTECT(R_fcall = lang2(func_musse, ScalarReal(*t)));
  pars = REAL(eval(R_fcall, rho_musse));
  do_derivs_musse(*neq / 2, pars, y, ydot);
  UNPROTECT(1);
}
