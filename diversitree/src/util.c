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
