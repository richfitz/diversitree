/* Multi-state BiSSE (MuSSE) compiled code */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>
#include "util.h"

/* For CVODES */
#include <nvector/nvector_serial.h>
#include <user_data.h>

/* This is the core function that actually evaluates the deriative */
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

  /* TODO: Replace with mult_mv() from util-matrix.c */
  /* mult_mv(k, Q, y, 1.0, ydot); */
  do_gemm2(Q, k, k, y, k, 2, ydot);
}

/* Plain MuSSE */
/* deSolve / LSODA */
static double *parms_musse;

void initmod_musse(void (* odeparms)(int *, double *)) {
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  parms_musse = REAL(get_deSolve_gparms());
} 

void derivs_musse(int *neq, double *t, double *y, double *ydot,
		   double *yout, int *ip) {
  do_derivs_musse(*neq / 2, parms_musse, y, ydot);
}

/* CVODES */
int derivs_musse_cvode(realtype t, N_Vector y, N_Vector ydot,
		       void *user_data) {
  const UserData *data = (UserData*) user_data;
  do_derivs_musse(data->neq/2,
		  data->p,
		  NV_DATA_S(y),
		  NV_DATA_S(ydot));
  return 0;
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

/* Time-dependent MuSSE */
/* The time dependent models are quite a bit different to the other
   models; rather than passing in a parameter vector, we pass in a (R)
   function of time.

   Each time step, evaluate this function to get the *actual*
   parameter vector, which is passed through to the underlying
   derivative calculation.

   The CVODES integrator requires that the model 'parameters' is a
   real vector, so we are going to pass around something trivial here
   instead.

   The actual function to be evaluated is tfunc_musse, and the required
   environment is trho_musse. */
static SEXP tfunc_musse;
static SEXP trho_musse;

void set_tfunc_musse_t(SEXP r_tfunc, SEXP r_trho) {
  if ( !isFunction(r_tfunc) )
    error("tfunc is not a function");
  if ( !isEnvironment(r_trho) )
    error("tenv is not an environment");

  tfunc_musse = r_tfunc;
  trho_musse  = r_trho;
}

void do_derivs_musse_t(int k, double t, double *y, double *ydot) {
  SEXP R_fcall, r_pars;
  double *pars;
  int i;
  const int np = k*(k+2), ncheck = 2*k;

  PROTECT(R_fcall = lang2(tfunc_musse, ScalarReal(t)));
  PROTECT(r_pars = eval(R_fcall, trho_musse));
  pars = REAL(r_pars);

  if ( LENGTH(r_pars) != np )
    error("Invalid parameter length");
  /* TODO: Only check for negative parameters in lambda and mu right
    now; the time varying Q matrix is still not done.  Once it is, we
    need to check the nondiagonal elements of the Q matrix */
  for ( i = 0; i < ncheck; i++ )
    if ( pars[i] < 0 )
      error("Illegal negative parameter at time %2.5f", t);

  do_derivs_musse(k, pars, y, ydot);

  UNPROTECT(2);
}

/* deSolve / LSODA */
void initmod_musse_t(void (* odeparms)(int *, double *)) {
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  SEXP obj = get_deSolve_gparms();
  set_tfunc_musse_t(VECTOR_ELT(obj, 0), VECTOR_ELT(obj, 1));
}

void derivs_musse_t(int *neq, double *t, double *y, double *ydot, 
		 double *yout, int *ip) {
  do_derivs_musse_t(*neq / 2, *t, y, ydot);
}

/* CVODES */
SEXP r_set_tfunc_musse_t(SEXP r_tfunc, SEXP r_trho) {
  set_tfunc_musse_t(r_tfunc, r_trho);
  return R_NilValue;
}

int derivs_musse_t_cvode(realtype t, N_Vector y, N_Vector ydot,
			 void *user_data) {
  do_derivs_musse_t(((UserData*) user_data)->neq/2,
		    t,
		    NV_DATA_S(y),
		    NV_DATA_S(ydot));
  return 0;
}
