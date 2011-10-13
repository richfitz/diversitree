/*
 * These are the BiSSE equations, implemented in c
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* For CVODES */
#include <nvector/nvector_serial.h>
#include <user_data.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_bisse(double *pars, double *y, double *ydot) {
  double E0 = y[0], E1 = y[1], D0 = y[2], D1 = y[3];
  double la0 = pars[0], la1 = pars[1], mu0 = pars[2], mu1=pars[3],
    q01 = pars[4], q10 = pars[5];

  ydot[0] = -(mu0 + q01 + la0) * E0 + la0 * E0 * E0 + mu0 + q01 * E1;
  ydot[1] = -(mu1 + q10 + la1) * E1 + la1 * E1 * E1 + mu1 + q10 * E0;
  ydot[2] = -(mu0 + q01 + la0) * D0 + 2 * la0 * E0 * D0 + q01 * D1;
  ydot[3] = -(mu1 + q10 + la1) * D1 + 2 * la1 * E1 * D1 + q10 * D0;
}

/* Plain BiSSE */
/* deSolve / LSODA */
static double parms_bisse[6];

void initmod_bisse(void (* odeparms)(int *, double *)) {
  int N = 6;
  odeparms(&N, parms_bisse);
}

void derivs_bisse(int *neq, double *t, double *y, double *ydot, 
		  double *yout, int *ip) {
  do_derivs_bisse(parms_bisse, y, ydot);
}

/* CVODES */
int derivs_bisse_cvode(realtype t, N_Vector y, N_Vector ydot,
		       void *user_data) {
  do_derivs_bisse(((UserData*) user_data)->p,
		  NV_DATA_S(y),
		  NV_DATA_S(ydot));
  return 0;
}

void initial_conditions_bisse(int neq, double *vars_l, double *vars_r,
			      double *pars, double t, 
			      double *vars_out) {
  vars_out[0] = vars_l[0];
  vars_out[1] = vars_l[1];
  vars_out[2] = vars_l[2] * vars_r[2] * pars[0];
  vars_out[3] = vars_l[3] * vars_r[3] * pars[1];
}

/* Not used, but this is the BiSSE Jacobian */
void do_jac_bisse(double *pars, double *y, double *pd, int n) {
  double E0 = y[0], E1 = y[1], D0 = y[2], D1 = y[3];
  double la0 = pars[0], la1 = pars[1], mu0 = pars[2], mu1=pars[3],
    q01 = pars[4], q10 = pars[5];

  pd[0] = -(mu0 + q01 + la0) + 2 * la0 * E0;
  pd[1] = q01;

  pd[n] = q10;
  pd[n+1] = -(mu1 + q10 + la1) + 2 * la1 * E1;

  pd[2*n] = 2 * D0 * la0;
  pd[2*n+2] = -(mu0 + q01 + la0) + 2 * la0 * E0;
  pd[2*n+3] = q01;

  pd[3*n+1] = 2 * D1 * la1;
  pd[3*n+2] = q01;
  pd[3*n+3] = -(mu1 + q10 + la1) + 2 * la1 * E1;
}

void jac_bisse(int *neq, double *t, double *y, int *ml, int *mu, 
	       double *pd, int *nrowpd, double *yout, int *ip) {
  do_jac_bisse(parms_bisse, y, pd, *nrowpd);
}

/* Time-dependant BiSSE */
/* The time dependent models are quite a bit different to the other
   models; rather than passing in a parameter vector, we pass in a (R)
   function of time.

   Each time step, evaluate this function to get the *actual*
   parameter vector, which is passed through to the underlying
   derivative calculation.

   The CVODES integrator requires that the model 'parameters' is a
   real vector, so we are going to pass around something trivial here
   instead.

   The actual function to be evaluated is tfunc_bisse, and the required
   environment is trho_bisse. */
static SEXP tfunc_bisse;
static SEXP trho_bisse;

void set_tfunc_bisse_t(SEXP r_tfunc, SEXP r_trho) {
  if ( !isFunction(r_tfunc) )
    error("tfunc is not a function");
  if ( !isEnvironment(r_trho) )
    error("tenv is not an environment");

  tfunc_bisse = r_tfunc;
  trho_bisse  = r_trho;
}

void do_derivs_bisse_t(double t, double *y, double *ydot) {
  SEXP R_fcall, r_pars;
  double *pars;
  int i;

  PROTECT(R_fcall = lang2(tfunc_bisse, ScalarReal(t)));
  PROTECT(r_pars = eval(R_fcall, trho_bisse));
  pars = REAL(r_pars);

  if ( LENGTH(r_pars) != 6 )
    error("Invalid parameter length");
  for ( i = 0; i < 6; i++ )
    if ( pars[i] < 0 )
      error("Illegal negative parameter at time %2.5f", t);

  do_derivs_bisse(pars, y, ydot);

  UNPROTECT(2);
}

/* deSolve / LSODA */
void initmod_bisse_t(void (* odeparms)(int *, double *)) {
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  SEXP obj = get_deSolve_gparms();
  set_tfunc_bisse_t(VECTOR_ELT(obj, 0), VECTOR_ELT(obj, 1));
}

void derivs_bisse_t(int *neq, double *t, double *y, double *ydot, 
		 double *yout, int *ip) {
  do_derivs_bisse_t(*t, y, ydot);
}

/* CVODES */
SEXP r_set_tfunc_bisse_t(SEXP r_tfunc, SEXP r_trho) {
  set_tfunc_bisse_t(r_tfunc, r_trho);
  return R_NilValue;
}

int derivs_bisse_t_cvode(realtype t, N_Vector y, N_Vector ydot,
		      void *user_data) {
  do_derivs_bisse_t(t,
		    NV_DATA_S(y),
		    NV_DATA_S(ydot));
  return 0;
}

/* Auxilliary (just compute E) */
void initmod_bisse_aux(void (* odeparms)(int *, double *)) {
  int N = 6;
  odeparms(&N, parms_bisse);
}

void do_derivs_bisse_aux(double *pars, double *y, double *ydot) {
  double E0 = y[0], E1 = y[1];
  const double la0 = pars[0], la1 = pars[1], 
    mu0 = pars[2], mu1=pars[3],
    q01 = pars[4], q10 = pars[5];

  ydot[0] = -(mu0 + q01 + la0) * E0 + la0 * E0 * E0 + mu0 + q01 * E1;
  ydot[1] = -(mu1 + q10 + la1) * E1 + la1 * E1 * E1 + mu1 + q10 * E0;
}

/* deSolve */
void derivs_bisse_aux(int *neq, double *t, double *y, double *ydot, 
		      double *yout, int *ip) {
  do_derivs_bisse_aux(parms_bisse, y, ydot);
}

/* CVODES */
int derivs_bisse_aux_cvode(realtype t, N_Vector y, N_Vector ydot,
			   void *user_data) {
  do_derivs_bisse_aux(((UserData*) user_data)->p,
		      NV_DATA_S(y),
		      NV_DATA_S(ydot));
  return 0;
}
