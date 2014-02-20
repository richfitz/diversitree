#include "config.h"

#ifdef HAVE_FFTW3_H

#include <R.h>
#include <Rdefines.h>
#include <complex.h>
#include <fftw3.h>
#include "rfftw.h"

rfftw_plan_real* make_rfftw_plan_real(int nd, int nx, int dir,
				      double *x, fftw_complex *y, 
				      int flags) {
  rfftw_plan_real *obj =
    (rfftw_plan_real*)calloc(1, sizeof(rfftw_plan_real));
  int ny = ((int)floor(nx/2)) + 1;
  int xstride, xdist, ystride, ydist;

  if ( dir == DIR_COLS ) {
    xstride = 1;
    ystride = 1;
    xdist   = nx;
    ydist   = ny;
  } else { /* if ( dir == DIR_ROWS ) { // assume ROWS silently */
    xstride = nd;
    ystride = nd;
    xdist   = 1;
    ydist   = 1;
  }

  obj->nd = nd;
  obj->nx = nx;
  obj->ny = ny;
  obj->x = x;
  obj->y = y;
  obj->dir = dir;

  /* Consider FFTW_PATIENT, or even FFTW_EXHAUSTIVE */
  obj->plan_f = fftw_plan_many_dft_r2c(1, &nx, nd, 
				       obj->x, NULL, xstride, xdist,
				       obj->y, NULL, ystride, ydist,
				       flags);
  obj->plan_b = fftw_plan_many_dft_c2r(1, &nx, nd, 
				       obj->y, NULL, ystride, ydist,
				       obj->x, NULL, xstride, xdist,
				       flags);
  return obj;
}

/* Next, the R interface */
static void rfftw_plan_real_finalize(SEXP extPtr);
SEXP r_make_rfftw_plan_real(SEXP r_nd, SEXP r_nx, SEXP r_dir) {
  const int flags = FFTW_PATIENT;
  SEXP extPtr;
  rfftw_plan_real *obj;  
  int nd, nx, ny, dir;
  double *x;
  fftw_complex *y;
  PROTECT(r_nd  = AS_INTEGER(r_nd));
  PROTECT(r_nx  = AS_INTEGER(r_nx));
  PROTECT(r_dir = AS_INTEGER(r_dir));

  nd = INTEGER(r_nd)[0];
  nx = INTEGER(r_nx)[0];
  dir = INTEGER(r_dir)[0];

  ny = ((int)floor(nx/2)) + 1;

  x = fftw_malloc(nd * nx * sizeof(double));	  
  y = fftw_malloc(nd * (ny+1) * sizeof(fftw_complex));

  obj = make_rfftw_plan_real(nd, nx, dir, x, y, flags);
  extPtr = R_MakeExternalPtr(obj, install("test_ptr"), R_NilValue);

  R_RegisterCFinalizer(extPtr, rfftw_plan_real_finalize);

  UNPROTECT(3);
  return extPtr;
}

static void rfftw_plan_real_finalize(SEXP extPtr) {
  rfftw_plan_real *obj =
    (rfftw_plan_real*)R_ExternalPtrAddr(extPtr);
  fftw_destroy_plan(obj->plan_f);
  fftw_destroy_plan(obj->plan_b);
  fftw_free(obj->x);
  fftw_free(obj->y);
  free(obj);
}

SEXP r_rfftw_forw(SEXP extPtr, SEXP r_x_in) {
  rfftw_plan_real *obj =
    (rfftw_plan_real*)R_ExternalPtrAddr(extPtr);
  int i, nxd, nyd;
  double  *x_in,  *x = obj->x;
  fftw_complex *y_out, *y = obj->y;
  SEXP ret;

  PROTECT(r_x_in = AS_NUMERIC(r_x_in));
  x_in = REAL(r_x_in);
  nxd = obj->nx * obj->nd;
  nyd = obj->ny * obj->nd;

  for ( i = 0; i < nxd; i++ )
    x[i] = x_in[i];

  fftw_execute(obj->plan_f);

  PROTECT(ret = allocVector(CPLXSXP, nyd));
  /* TODO: There is a change that I should be using Rcomplex here */
  y_out = (fftw_complex*)COMPLEX(ret);
  for ( i = 0; i < nyd; i++ )
    y_out[i] = y[i];

  UNPROTECT(2);
  return ret;
}

SEXP r_rfftw_back(SEXP extPtr, SEXP r_y_in) {
  rfftw_plan_real *obj = 
    (rfftw_plan_real*)R_ExternalPtrAddr(extPtr);
  int i, nxd, nyd;
  fftw_complex *y_in,  *y = obj->y;
  double  *x_out, *x = obj->x;
  SEXP ret;
  
  PROTECT(r_y_in = AS_COMPLEX(r_y_in));
  /* TODO: There is a change that I should be using Rcomplex here */
  y_in = (fftw_complex*)COMPLEX(r_y_in);
  nxd = obj->nx * obj->nd;
  nyd = obj->ny * obj->nd;

  for ( i = 0; i < nyd; i++ )
    y[i] = y_in[i];

  fftw_execute(obj->plan_b);

  PROTECT(ret = allocVector(REALSXP, nxd));
  x_out = REAL(ret);
  for ( i = 0; i < nxd; i++ )
    x_out[i] = x[i];

  UNPROTECT(2);
  return ret;
}

/* Simple help for getting and setting wisdom */
SEXP r_get_wisdom() {
  char *wisdom = fftw_export_wisdom_to_string();
  SEXP ret;
  PROTECT(ret = allocVector(STRSXP, 1));
  SET_STRING_ELT(ret, 0, mkChar(wisdom));
  UNPROTECT(1);
  return ret;
}

SEXP r_set_wisdom(SEXP r_wisdom) {
  const char *wisdom = CHAR(STRING_ELT(r_wisdom, 0));
  SEXP ret;
  PROTECT(ret = allocVector(LGLSXP, 1));
  INTEGER(ret)[0] = fftw_import_wisdom_from_string(wisdom);
  UNPROTECT(1);
  return ret;
}

#endif
