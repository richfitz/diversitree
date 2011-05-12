#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

/* Simplified matrix multiplication, assuming straightforward sizes
   and zeroing the input.  GEMM does:
     Z = alpha X Y + beta Z

   With 
     beta = 0, this gives a fresh calculation of X Y, and with 
     beta = 1, this updates Z to give Z = Z + X Y
*/
void do_gemm(double *x, int nrx, int ncx,
             double *y, int nry, int ncy,
             double *z) {
  const char *trans = "N";
  double alpha = 1.0, beta = 0.0;
  F77_CALL(dgemm)(trans, trans, &nrx, &ncy, &ncx, &alpha,
                  x, &nrx, y, &nry, &beta, z, &nrx);
}

void do_gemm2(double *x, int nrx, int ncx,
	      double *y, int nry, int ncy,
	      double *z) {
  const char *trans = "N";
  double alpha = 1.0, beta = 1.0;
  F77_CALL(dgemm)(trans, trans, &nrx, &ncy, &ncx, &alpha,
                  x, &nrx, y, &nry, &beta, z, &nrx);
}

void do_gemm3(double *x, int nrx, int ncx,
	      double *y, int nry, int ncy,
	      double *z, double beta) {
  const char *trans = "N";
  double alpha = 1.0;
  F77_CALL(dgemm)(trans, trans, &nrx, &ncy, &ncx, &alpha,
                  x, &nrx, y, &nry, &beta, z, &nrx);
}

void r_gemm(double *x, int *nrx, int *ncx,
            double *y, int *nry, int *ncy,
            double *z) {
  do_gemm(x, *nrx, *ncx, y, *nry, *ncy, z);
}

void r_gemm2(double *x, int *nrx, int *ncx,
	     double *y, int *nry, int *ncy,
	     double *z) {
  do_gemm2(x, *nrx, *ncx, y, *nry, *ncy, z);
}

SEXP matrix_to_list(SEXP r_m) {
  SEXP ret, tmp;
  int i, j, k, nr = nrows(r_m), nc = ncols(r_m);
  double *in, *out;

  in = REAL(r_m);

  PROTECT(ret = allocVector(VECSXP, nr));

  for ( i = 0; i < nr; i++ ) {
    /*
      I believe that I don't have to protect agressively here;
      otherwise something like below would be needed.

      PROTECT(tmp = allocVector(REALSXP, nc));
      SET_VECTOR_ELT(ret, i, tmp);
      UNPROTECT(1);
      out = REAL(tmp);

      another option, which definitely does not need garbage
      collection, is:

      SET_VECTOR_ELT(ret, i, allocVector(REALSXP, nc));
      out = REAL(VECTOR_ELT(ret, i));

      which falls somewhere between the two approaches in speed.
      
      However, I've run this under gctorture, and it seems not to
      crash, which is a good sign.
    */
    SET_VECTOR_ELT(ret, i, tmp = allocVector(REALSXP, nc));
    out = REAL(tmp);

    for ( j = 0, k = i; j < nc; j++, k+= nr )
      out[j] = in[k];
  }

  UNPROTECT(1);
  return ret;
}
