#include <R.h>
#include <Rinternals.h>
#include "util-splines.h"

/* This is adapted directly from the R sources:
   
   ${R_SRC}/src/library/stats/src/splines.c

   As this is not API, we can't use the functions directly from C.
   However, I need a simple set of spline functions for use when
   computing derivatives without calling R functions.

   This strips down the file [describe].

   The original copyright notice follows below.  I claim no authorship
   over this piece of work.
*/

/* Note that 'n' must be a single integer */
void RSRC_fmm_spline(int n, double *x, double *y, 
		     double *b, double *c, double *d) {
  int nm1, i;
  double t;

  /* Adjustment for 1-based arrays */

  x--; y--; b--; c--; d--;

  if(n < 2) {
    error("Too few points");
  }

  if(n < 3) {
    t = (y[2] - y[1]);
    b[1] = t / (x[2]-x[1]);
    b[2] = b[1];
    c[1] = c[2] = d[1] = d[2] = 0.0;
    return;
  }

  nm1 = n - 1;

  /* Set up tridiagonal system */
  /* b = diagonal, d = offdiagonal, c = right hand side */

  d[1] = x[2] - x[1];
  c[2] = (y[2] - y[1])/d[1];/* = +/- Inf	for x[1]=x[2] -- problem? */
  for(i=2 ; i<n ; i++) {
    d[i] = x[i+1] - x[i];
    b[i] = 2.0 * (d[i-1] + d[i]);
    c[i+1] = (y[i+1] - y[i])/d[i];
    c[i] = c[i+1] - c[i];
  }

  /* End conditions. */
  /* Third derivatives at x[0] and x[n-1] obtained */
  /* from divided differences */

  b[1] = -d[1];
  b[n] = -d[nm1];
  c[1] = c[n] = 0.0;
  if(n > 3) {
    c[1] = c[3]/(x[4]-x[2]) - c[2]/(x[3]-x[1]);
    c[n] = c[nm1]/(x[n] - x[n-2]) - c[n-2]/(x[nm1]-x[n-3]);
    c[1] = c[1]*d[1]*d[1]/(x[4]-x[1]);
    c[n] = -c[n]*d[nm1]*d[nm1]/(x[n]-x[n-3]);
  }

  /* Gaussian elimination */

  for(i=2 ; i<=n ; i++) {
    t = d[i-1]/b[i-1];
    b[i] = b[i] - t*d[i-1];
    c[i] = c[i] - t*c[i-1];
  }

  /* Backward substitution */

  c[n] = c[n]/b[n];
  for(i=nm1 ; i>=1 ; i--)
    c[i] = (c[i]-d[i]*c[i+1])/b[i];

  /* c[i] is now the sigma[i-1] of the text */
  /* Compute polynomial coefficients */

  b[n] = (y[n] - y[n-1])/d[n-1] + d[n-1]*(c[n-1]+ 2.0*c[n]);
  for(i=1 ; i<=nm1 ; i++) {
    b[i] = (y[i+1]-y[i])/d[i] - d[i]*(c[i+1]+2.0*c[i]);
    d[i] = (c[i+1]-c[i])/d[i];
    c[i] = 3.0*c[i];
  }
  c[n] = 3.0*c[n];
  d[n] = d[nm1];
  return;
}

void RSRC_fmm_spline_eval(int nu, double *u, double *v,
			  int n,  double *x, double *y, 
			  double *b, double *c, double *d) {
  /* Evaluate  v[l] := spline(u[l], ...),	    l = 1,..,nu, i.e. 0:(nu-1)
   * Nodes x[i], coef (y[i]; b[i],c[i],d[i]); i = 1,..,n , i.e. 0:(n-1)
   */
  const int n_1 = n - 1;
  int i, j, k, l;
  double ul, dx, tmp;

  for(l = 0; l < nu; l++)
    v[l] = u[l];

  i = 0;
  for(l = 0; l < nu; l++) {
    ul = v[l];
    if(ul < x[i] || (i < n_1 && x[i+1] < ul)) {
      /* reset i  such that  x[i] <= ul <= x[i+1] : */
      i = 0;
      j = n;
      do {
	k = (i+j)/2;
	if(ul < x[k]) j = k;
	else i = k;
      }
      while(j > i+1);
    }
    dx = ul - x[i];

    tmp = d[i];

    v[l] = y[i] + dx*(b[i] + dx*(c[i] + dx*tmp));
  }
}


/* Below here, written by RGF */

dt_spline* make_dt_spline(int nx, double *x, double *y, int deriv) {
  double *b, *c, *d;
  int i;
  dt_spline *obj = (dt_spline *)Calloc(1, dt_spline);
  obj->nx = nx;

  obj->x     = (double*) Calloc(nx, double);
  obj->y     = (double*) Calloc(nx, double);
  obj->b = b = (double*) Calloc(nx, double);
  obj->c = c = (double*) Calloc(nx, double);
  obj->d = d = (double*) Calloc(nx, double);

  memcpy(obj->x, x, nx * sizeof(double));
  memcpy(obj->y, y, nx * sizeof(double));

  RSRC_fmm_spline(nx, obj->x, obj->y, b, c, d);
  
  if ( deriv > 0 ) {
    for ( i = 0; i < nx; i++ ) {
      obj->y[i] =   b[i]; /* Alter y, use b */
      b[i]      = 2*c[i]; /* Alter b, use c */
      c[i]      = 3*d[i]; /* Alter c, use d */
      d[i]      =   0;    /* Alter d, use 0 */
    }
  }

  return obj;
}

void cleanup_dt_spline(dt_spline *obj) {
  Free(obj->x);
  Free(obj->y);
  Free(obj->b);
  Free(obj->c);
  Free(obj->d);

  Free(obj);
}

static void dt_spline_finalize(SEXP extPtr);
SEXP r_make_dt_spline(SEXP x, SEXP y, SEXP r_deriv) {
  SEXP extPtr;
  dt_spline *obj = 
    make_dt_spline(LENGTH(x), REAL(x), REAL(y), INTEGER(r_deriv)[0]);
  extPtr = R_MakeExternalPtr(obj, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(extPtr, dt_spline_finalize);
  return extPtr;
}

static void dt_spline_finalize(SEXP extPtr) {
  cleanup_dt_spline((dt_spline*)R_ExternalPtrAddr(extPtr));
}

/* A couple of different evaluation functions */
double dt_spline_eval1(dt_spline *obj, double u) {
  double ret;
  RSRC_fmm_spline_eval(1, &u, &ret, 
		       obj->nx, obj->x, obj->y, 
		       obj->b, obj->c, obj->d);
  return ret;
}

void dt_spline_eval(dt_spline *obj, double *u, int nu, double *ret) {
  RSRC_fmm_spline_eval(nu, u, ret, 
		       obj->nx, obj->x, obj->y, 
		       obj->b, obj->c, obj->d);
}

/* From R, evaluate points on this spline */
SEXP r_dt_spline_eval(SEXP extPtr, SEXP u) {
  int nu = LENGTH(u);
  dt_spline *obj = (dt_spline*)R_ExternalPtrAddr(extPtr);
  SEXP ret;

  PROTECT(ret = allocVector(REALSXP, nu));
  dt_spline_eval(obj, REAL(u), nu, REAL(ret));
  UNPROTECT(1);
  return ret;
}
