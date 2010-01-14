/* Mkn - Markov k-state n-parameter character model */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>

static double *parms_mkn;
void initmod_mkn(void (* odeparms)(int *, double *)) {
  /* TODO: I should check here about the lengths of parameters, but I
     won't bother; it is not clear how best to do this, anyway.  Most
     of the checking should be done in the R end of things.  Because
     this is going to get an R object, it should be much easier to the
     checking there. */
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  parms_mkn = REAL(get_deSolve_gparms());
} 

/* Simplified matrix multiplication, assuming straightforward sizes
   and zeroing the input.  GEMM does:
     Z = alpha X Y + beta Z
*/
void do_gemm(double *x, int nrx, int ncx,
             double *y, int nry, int ncy,
             double *z) {
  const char *trans = "N";
  double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)(trans, trans, &nrx, &ncy, &ncx, &one,
                  x, &nrx, y, &nry, &zero, z, &nrx);
}

void derivs_mkn(int *neq, double *t, double *y, double *ydot,
		double *yout, int *ip) {
  const int k = *neq;
  do_gemm(parms_mkn, k, k, y, k, 1, ydot);
}

void r_gemm(double *x, int *nrx, int *ncx,
            double *y, int *nry, int *ncy,
            double *z) {
  do_gemm(x, *nrx, *ncx, y, *nry, *ncy, z);
}
