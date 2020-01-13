/* Utilities for working with matrix multiplication */

#define USE_FC_LEN_T
#include <Rconfig.h>
#include <R_ext/BLAS.h>
#ifndef FCONE
# define FCONE
#endif

#include <R.h>
#include "util-matrix.h"

/* The routines are 
     [z]mult_ab[c]
   where
     z (if present) indicates all types are Rcomplex
     a, b, and (if present) c are
       m: dense matrix
       d: vector indicating elements of a diagonal matrix
       v: vector
   All matrices are assumed to be square.

   Only a subset of possible routines are included below.

   Only mult_mv() includes DGEMM's "beta" argument for incrementing a
   product.
*/

/* Compute 
     out = A * v + beta * out 
   for a square n x n matrix A and a n-long vector v.  'beta' will
   typically be 0.0 or 1.0.
*/
void mult_mv(int n, double *A, double *v, double beta, double *out) {
  const char *trans = "N";
  int one = 1;
  double alpha = 1.0;
  F77_CALL(dgemm)(trans, trans, &n, &one, &n, &alpha,
                  A, &n, v, &n, &beta, out, &n FCONE FCONE);
}

/* Compute out = A * B for two n x n matrices */
void mult_mm(int n, double *A, double *B, double *out) {
  const char *trans = "N";
  double alpha = 1.0, beta = 0.0;
  F77_CALL(dgemm)(trans, trans, &n, &n, &n, &alpha,
                  A, &n, B, &n, &beta, out, &n FCONE FCONE);
}

/* This computes 
   A * d
   where A is an arbitrary (nonsymmetric) n*n dense matrix and d is a
   n-long vector representing the diagonal elements of a matrix D */
void mult_md(int n, double *A, double *d, double *out) {
  int i, j, k=0;
  for ( i = 0; i < n; i++ )
    for ( j = 0; j < n; j++, k++ )
      out[k] = d[i] * A[k];
}

/* Compute out = A * B * C for three n x n matrices.  The space 'wrk'
   must be at least n*n long.
*/
void mult_mmm(int n, double *A, double *B, double *C, 
	      double *out, double *wrk) {
  mult_mm(n, A, B, wrk);
  mult_mm(n, wrk, C, out);
}

/* Compute out = A * d * C for two n x n matrices (A and C), and a
   diagonal matrix whose n elements are given by the vector d.  The
   space 'wrk' must be at least n*n long.
*/
void mult_mdm(int n, double *A, double *d, double *C, double *out,
	      double *wrk) {
  mult_md(n, A, d, wrk);
  mult_mm(n, wrk, C, out);
}

/* R interfaces */
void r_mult_mv(int *n, double *A, double *v, double *beta, double *out) {
  mult_mv(*n, A, v, *beta, out);
}

void r_mult_mm(int *n, double *A, double *B, double *out) {
  mult_mm(*n, A, B, out);
}

void r_mult_md(int *n, double *A, double *d, double *out) {
  mult_md(*n, A, d, out);
}

void r_mult_mmm(int *n, double *A, double *B, double *C, double *out,
		double *wrk) {
  mult_mmm(*n, A, B, C, out, wrk);
}

void r_mult_mdm(int *n, double *A, double *d, double *C, double *out,
		double *wrk) {
  mult_mdm(*n, A, d, C, out, wrk);
}
