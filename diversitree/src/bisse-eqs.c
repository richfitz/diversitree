/*
 * These are the BiSSE equations, implemented in c
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static SEXP func_bisse;
static SEXP rho_bisse;
static double parms_bisse[6];

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
void initmod_bisse(void (* odeparms)(int *, double *)) {
  int N = 6;
  odeparms(&N, parms_bisse);
}

void derivs_bisse(int *neq, double *t, double *y, double *ydot, 
		  double *yout, int *ip) {
  do_derivs_bisse(parms_bisse, y, ydot);
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
void initmod_bisse_t(void (* odeparms)(int *, double *)) {
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  SEXP obj = get_deSolve_gparms();

  func_bisse  = VECTOR_ELT(obj, 0);
  rho_bisse   = VECTOR_ELT(obj, 1);

  if ( !isFunction(func_bisse) )
    error("First element must be a function");
  if ( !isEnvironment(rho_bisse) )
    error("Second element must be an environment");
}

void derivs_bisse_t(int *neq, double *t, double *y, double *ydot, 
		    double *yout, int *ip) {
  SEXP R_fcall, r_pars;
  double *pars;

  PROTECT(R_fcall = lang2(func_bisse, ScalarReal(*t)));
  r_pars = eval(R_fcall, rho_bisse);
  pars = REAL(r_pars);

  if ( LENGTH(r_pars) != 6 )
    error("Invalid parameter length");
  /* TODO:
  for ( i = 0; i < 6; i++ )
    if ( variable[i] && pars[i] < 0 )
      error("Illegal negative parameter");
  */

  do_derivs_bisse(pars, y, ydot);

  UNPROTECT(1);
}
