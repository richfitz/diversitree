/* Simplified routines for a few linear algebra things */
/* TODO: A function that takes the returned elements of la_eigen() or
   la_eigen_oneshot() and converts them to Rcomplex types */
#include <R.h>
#include <R_ext/Lapack.h>

/* 1: invert a matrix */
/* This is destructive on A, so a copy needs to be made */
/* ipiv needs to be at least 'n' long, and 'wrk' needs to be at least
   n*n long. */
void la_invert(int n, double *A, int *ipiv, double *wrk) {
  int lwork = n * n;
  int info;
  F77_NAME(dgetrf)(&n, &n, A, &n, ipiv, &info);
  F77_NAME(dgetri)(&n, A, &n, ipiv, wrk, &lwork, &info);
}

/* 2: Compute the eigenvalues and right eigenvectors of a matrix */
/*
  In what follows, the work vectors "wr" and "wi" must be at least
  length N, "vl" and "vr" are length N^2 (these will store the
  eigenvectors).  But if we don't compute the left eigenvectors then
  vl can be NULL, I believe.
*/

#define VECTORS "Vectors"
#define NOVECTORS "No Vectors"
/* a: determine the size of the work array */
/*   n = dimension of A
     A = square matrix
     values = eigenvalues
     values_i = imaginary part of the eigenvalues
     vectors = right eigenvectors
*/
int la_eigen_query(int n, double *A, double *values, double *values_i,
		   double *vectors) {
  double wkopt;
  int lwork = -1, info;
  double *vl = NULL;
  F77_NAME(dgeev)(NOVECTORS, VECTORS, &n, A, &n, values, values_i, vl,
		  &n, vectors, &n, &wkopt, &lwork, &info);
  lwork = (int)wkopt;
  if ( info != 0 )
    error("Error determining space for eigenproblem");
  return lwork;
}

/* b: Actually compute the problem*/
/* This is destructive on A, so a copy needs to be made */
int la_eigen(int n, double *A, double *values, double *values_i, 
	     double *vectors, int lwork, double *work) {
  int info;
  double *vl = NULL;
  F77_NAME(dgeev)(NOVECTORS, VECTORS, &n, A, &n, values, values_i, vl, &n, vectors,
		  &n, work, &lwork, &info);
  return info;
}

/* c: one shot version that does dynamic memory allocation */
int la_eigen_oneshot(int n, double *A, double *values, 
		     double *values_i, double *vectors) {
  int lwork, info;
  double *work;
  
  lwork = la_eigen_query(n, A, values, values_i, vectors);
  work = (double *)R_alloc(lwork, sizeof(double));
  info = la_eigen(n, A, values, values_i, vectors, lwork, work);
  if ( info != 0 )
    error("Error computing eigensystem");
  return info;
}

/* R intefaces */
void r_la_invert(int *n, double *A, int *ipiv, double *wrk) {
  la_invert(*n, A, ipiv, wrk);
}

void r_la_eigen_query(int *n, double *A, double *values, 
		      double *values_i, double *vectors) {
  la_eigen_query(*n, A, values, values_i, vectors);
}

void r_la_eigen(int *n, double *A, double *values, double *values_i,
		double *vectors, int *lwork, double *work) {
  la_eigen(*n, A, values, values_i, vectors, *lwork, work);
}

void r_la_eigen_oneshot(int *n, double *A, double *values,
			double *values_i, double *vectors) {
  la_eigen_oneshot(*n, A, values, values_i, vectors);
}

