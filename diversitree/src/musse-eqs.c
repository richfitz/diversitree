/* Multi-state BiSSE compiled code */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>

void do_gemm2(double *x, int nrx, int ncx,
	      double *y, int nry, int ncy,
	      double *z) {
  const char *trans = "N";
  double alpha = 1.0, beta = 1.0;
  F77_CALL(dgemm)(trans, trans, &nrx, &ncy, &ncx, &alpha,
                  x, &nrx, y, &nry, &beta, z, &nrx);
}

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

/* This gives a skeleton for the final version */
void derivs_musse(int *neq, double *t, double *y, double *ydot,
		   double *yout, int *ip) {
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
}

void r_gemm2(double *x, int *nrx, int *ncx,
	     double *y, int *nry, int *ncy,
	     double *z) {
  do_gemm2(x, *nrx, *ncx, y, *nry, *ncy, z);
}

void do_gemm3(double *x, int nrx, int ncx,
	      double *y, int nry, int ncy,
	      double *z, double *beta) {
  const char *trans = "N";
  double alpha = 1.0;
  F77_CALL(dgemm)(trans, trans, /* trans = "N" - dont transpose */
		  &nrx,   /* M */
		  &ncy,   /* N */
		  &ncx,   /* K */
		  &alpha, /* alpha */
                  x,      /* A */
		  &nrx,   /* LDA */
		  y,      /* B */
		  &nry,   /* LDB */
		  beta,  /* beta */
		  z,      /* C */
		  &nrx);  /* LDC */
}

void r_gemm3(double *x, int *nrx, int *ncx,
	     double *y, int *nry, int *ncy,
	     double *z, double *beta) {
  do_gemm3(x, *nrx, *ncx, y, *nry, *ncy, z, beta);
}



/* This gives a skeleton for the final version */
void derivs_musse2(double *pars,
		   int *neq, double *t, double *y, double *ydot,
		   double *yout, int *ip) {
  const int k = *neq / 2;
  double *e = y, *d = y + k;
  double *dEdt = ydot, *dDdt = ydot + k;
  double *lambda = pars, *mu = pars + k, *Q = pars + 2*k;
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
}
