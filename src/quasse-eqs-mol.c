/* Code for doing MOL style integration.  The main headache here is
   that the amount of information that must be copied across is fairly
   high */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>

/* Here, 'diffusion' must be is transformed before being passed in: it
   incorporates the "diff/2" term that appears in the QuaSSE equation
   and the 1/dx^2 term that appears in the FD approximation to the
   second derivative.
*/

static double *parms_quasse;
void initmod_quasse_mol(void (* odeparms)(int *, double *)) {
  /* TODO: I should check here about the lengths of parameters, but I
     won't bother; it is not clear how best to do this, anyway.  Most
     of the checking should be done in the R end of things.  Because
     this is going to get an R object, it should be much easier to the
     checking there. */
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  parms_quasse = REAL(get_deSolve_gparms());
} 

/* Note that drift is not currently supported */
void derivs_quasse_mol(int *neq, double *t, double *y, double *ydot, 
		       double *yout, int *ip) {
  const int nx = *neq / 2;

  double *E = y, *D = y + nx;
  double *E0 = E-1, *E1 = E, *E2 = E+1;
  double *D0 = D-1, *D1 = D, *D2 = D+1;

  double *dEdt = ydot, *dDdt = ydot + nx;

  double *lambda = parms_quasse, *mu = parms_quasse + nx,
    *lambda_mu = parms_quasse + 2*nx,
    /*drift = parms_quasse[3*nx], */ diffusion = parms_quasse[3*nx+1];

  int j;

  /* Diffusion term: E */
  dEdt[0]    = (E1[1] - E1[0]) * diffusion;
  for ( j = 1; j < (nx-1); j++ )
    dEdt[j] = (E0[j] - 2*E1[j] + E2[j]) * diffusion;
  dEdt[nx-1] = (E1[nx-2] - E1[nx-1]) * diffusion;

  /* D */
  dDdt[0]    = (D1[1] - D1[0]) * diffusion;
  for ( j = 1; j < (nx-1); j++ )
    dDdt[j] = (D0[j] - 2*D1[j] + D2[j]) * diffusion;
  dDdt[nx-1] = (D1[nx-2] - D1[nx-1]) * diffusion;

  /* Source terms */
  for ( j = 0; j < nx; j++ )
    dEdt[j] += mu[j] - lambda_mu[j]*E[j] + lambda[j]*E[j]*E[j];

  for ( j = 0; j < nx; j++ )
    dDdt[j] += -lambda_mu[j]*D[j] + 2*lambda[j]*D[j]*E[j];
}

