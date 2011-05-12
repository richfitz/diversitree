/*
  Constant rate Birth-death functions
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* For CVODES */
#include <nvector/nvector_serial.h>
#include <user_data.h>

/* This is the core function that actually evaluates the deriative */
void do_derivs_bd(double *pars, double *y, double *ydot) {
  double E = y[0], D = y[1];
  double lambda = pars[0], mu = pars[1];

  ydot[0] = mu - (mu + lambda)*E +   lambda*E*E;
  ydot[1] =    - (mu + lambda)*D + 2*lambda*D*E;
}

/* Plain BD */
/* deSolve / LSODA */
static double parms_bd[2];

void initmod_bd(void (* odeparms)(int *, double *)) {
  int N = 2;
  odeparms(&N, parms_bd);
}

void derivs_bd(int *neq, double *t, double *y, double *ydot, 
	       double *yout, int *ip) {
  do_derivs_bd(parms_bd, y, ydot);
}

/* CVODES */
int derivs_bd_cvode(realtype t, N_Vector y, N_Vector ydot,
		    void *user_data) {
  do_derivs_bd(((UserData*) user_data)->p,
	       NV_DATA_S(y),
	       NV_DATA_S(ydot));
  return 0;
}

void initial_conditions_bd(int neq, double *vars_l, double *vars_r,
			   double *pars, double t, 
			   double *vars_out) {
  vars_out[0] = vars_l[0];
  vars_out[1] = vars_l[1] * vars_r[1] * pars[0];
}

/* Time-dependant BD */
/* The time dependent models are quite a bit different to the other
   models; rather than passing in a parameter vector, we pass in a (R)
   function of time.

   Each time step, evaluate this function to get the *actual*
   parameter vector, which is passed through to the underlying
   derivative calculation.

   The CVODES integrator requires that the model 'parameters' is a
   real vector, so we are going to pass around something trivial here
   instead.

   The actual function to be evaluated is tfunc_bd, and the required
   environment is trho_bd. */
static SEXP tfunc_bd;
static SEXP trho_bd;

void set_tfunc_bd_t(SEXP r_tfunc, SEXP r_trho) {
  if ( !isFunction(r_tfunc) )
    error("tfunc is not a function");
  if ( !isEnvironment(r_trho) )
    error("tenv is not an environment");

  tfunc_bd = r_tfunc;
  trho_bd  = r_trho;
}

void do_derivs_bd_t(double t, double *y, double *ydot) {
  SEXP R_fcall, r_pars;
  double *pars;

  PROTECT(R_fcall = lang2(tfunc_bd, ScalarReal(t)));
  PROTECT(r_pars = eval(R_fcall, trho_bd));
  pars = REAL(r_pars);

  if ( LENGTH(r_pars) != 2 )
    error("Invalid parameter length");
  if ( pars[0] < 0 || pars[1] < 0 )
    error("Illegal negative parameter at time %2.5f", t);

  do_derivs_bd(pars, y, ydot);

  UNPROTECT(2);
}

/* deSolve / LSODA */
void initmod_bd_t(void (* odeparms)(int *, double *)) {
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  SEXP obj = get_deSolve_gparms();
  set_tfunc_bd_t(VECTOR_ELT(obj, 0), VECTOR_ELT(obj, 1));
}

void derivs_bd_t(int *neq, double *t, double *y, double *ydot, 
		 double *yout, int *ip) {
  do_derivs_bd_t(*t, y, ydot);
}

/* CVODES */
SEXP r_set_tfunc_bd_t(SEXP r_tfunc, SEXP r_trho) {
  set_tfunc_bd_t(r_tfunc, r_trho);
  return R_NilValue;
}

int derivs_bd_t_cvode(realtype t, N_Vector y, N_Vector ydot,
		      void *user_data) {
  do_derivs_bd_t(t,
		 NV_DATA_S(y),
		 NV_DATA_S(ydot));
  return 0;
}
