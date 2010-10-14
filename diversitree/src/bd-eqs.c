/*
  Constant rate Birth-death functions
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


static SEXP func_bd;
static SEXP rho_bd;
static double parms_bd[2];

void do_derivs_bd(double *pars, double *y, double *ydot) {
  double E = y[0], D = y[1];
  double lambda = pars[0], mu = pars[1];

  ydot[0] = mu - (mu + lambda)*E +   lambda*E*E;
  ydot[1] =    - (mu + lambda)*D + 2*lambda*D*E;
}

/* Plain BD */
void initmod_bd(void (* odeparms)(int *, double *)) {
  int N = 2;
  odeparms(&N, parms_bd);
}

void derivs_bd(int *neq, double *t, double *y, double *ydot, 
	       double *yout, int *ip) {
  do_derivs_bd(parms_bd, y, ydot);
}

/* This is identical (copy/paste) to initmod_bisse_t(), but replacing
   func_bisse with func_bd (and same with rho).  Doing this "properly"
   in C probably requires use of macros */
void initmod_bd_t(void (* odeparms)(int *, double *)) {
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  SEXP obj = get_deSolve_gparms();

  func_bd  = VECTOR_ELT(obj, 0);
  rho_bd   = VECTOR_ELT(obj, 1);

  if ( !isFunction(func_bd) )
    error("First element must be a function");
  if ( !isEnvironment(rho_bd) )
    error("Second element must be an environment");
}

/* This is also almost identical to the BiSSE analogue, but replacing
   the parameter length check with 2 */
void derivs_bd_t(int *neq, double *t, double *y, double *ydot, 
		    double *yout, int *ip) {
  SEXP R_fcall, r_pars;
  double *pars;

  PROTECT(R_fcall = lang2(func_bd, ScalarReal(*t)));
  r_pars = eval(R_fcall, rho_bd);
  pars = REAL(r_pars);

  if ( LENGTH(r_pars) != 2 )
    error("Invalid parameter length");
  /* TODO:
  for ( i = 0; i < 2; i++ )
    if ( variable[i] && pars[i] < 0 )
      error("Illegal negative parameter");
  */

  do_derivs_bd(pars, y, ydot);

  UNPROTECT(1);
}
