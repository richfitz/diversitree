#include "config.h"
#ifdef HAVE_FFTW3_H

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h> /* why? */
#include <Rmath.h> /* for dnorm() */
#include <complex.h>
#include <fftw3.h>
#include "rfftw.h"
#include "quasse-eqs-fftC.h"

static void quasse_fft_finalize(SEXP extPtr);

SEXP r_make_quasse_fft(SEXP r_nx, SEXP r_dx, SEXP r_nd, SEXP r_flags) {
  quasse_fft *obj;
  SEXP extPtr;
  int nx = INTEGER(r_nx)[0];
  double dx = REAL(r_dx)[0];
  int n_fft = LENGTH(r_nd);
  int i;
  int flags;
  int *nd = (int*)calloc(n_fft, sizeof(int));
  for ( i = 0; i < n_fft; i++ )
    nd[i] = INTEGER(r_nd)[i];
  
  /* Simple interface to FFTW's flags */
  if ( INTEGER(r_flags)[0] == -1 )
    flags = FFTW_ESTIMATE;
  else if ( INTEGER(r_flags)[0] == 1 )
    flags = FFTW_PATIENT;
  else if ( INTEGER(r_flags)[0] == 2 )
    flags = FFTW_EXHAUSTIVE;
  else
    flags = FFTW_MEASURE;

  obj = make_quasse_fft(n_fft, nx, dx, nd, flags);

  extPtr = R_MakeExternalPtr(obj, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(extPtr, quasse_fft_finalize);

  return extPtr;
}

SEXP r_set_x(SEXP extPtr, SEXP x) {
  quasse_fft *obj = (quasse_fft*)R_ExternalPtrAddr(extPtr);
  int nd = LENGTH(x) / obj->nx;
  qf_copy_x(obj, REAL(x), nd, 1);
  return R_NilValue;
}

SEXP r_get_x(SEXP extPtr, SEXP r_nd) {
  quasse_fft *obj = (quasse_fft*)R_ExternalPtrAddr(extPtr);
  SEXP x;
  int nd = INTEGER(r_nd)[0];
  PROTECT(x = allocMatrix(REALSXP, obj->nx, nd));
  qf_copy_x(obj, REAL(x), nd, 0);
  UNPROTECT(1);
  return x;
}

/* Only for debugging */
SEXP r_propagate_t(SEXP extPtr, SEXP vars, SEXP lambda, SEXP mu, SEXP dt) {
  SEXP ret;
  quasse_fft *obj = (quasse_fft*)R_ExternalPtrAddr(extPtr);
  int i, ndat=LENGTH(lambda);
  double c_dt = REAL(dt)[0];
  int idx=-1, nd = LENGTH(vars) / obj->nx;

  idx = lookup(nd, obj->nd, obj->n_fft);
  if ( idx < 0 )
    error("Failed to find nd = %d\n", nd);

  qf_copy_x(obj, REAL(vars), nd, 1);

  obj->lambda = REAL(lambda);
  obj->mu = REAL(mu);
  obj->ndat = ndat;

  for ( i = 0; i < ndat; i++ )
    obj->z[i] = exp(c_dt * (obj->lambda[i] - obj->mu[i]));

  propagate_t(obj, idx);

  obj->lambda = NULL;
  obj->mu = NULL;

  PROTECT(ret = allocMatrix(REALSXP, obj->nx, nd));
  qf_copy_x(obj, REAL(ret), nd, 0);
  UNPROTECT(1);
  return ret;
}

SEXP r_propagate_x(SEXP extPtr, SEXP vars, SEXP drift, SEXP diffusion,
		   SEXP dt, SEXP padding) {
  SEXP ret;
  quasse_fft *obj = (quasse_fft*)R_ExternalPtrAddr(extPtr);
  int nkl = INTEGER(padding)[0], nkr = INTEGER(padding)[1];
  int idx=-1, nd=LENGTH(vars) / obj->nx;

  idx = lookup(nd, obj->nd, obj->n_fft);
  if ( idx < 0 )
    error("Failed to find nd = %d\n", nd);

  qf_copy_x(obj, REAL(vars), nd, 1);
  qf_setup_kern(obj, REAL(drift)[0], REAL(diffusion)[0], REAL(dt)[0],
		nkl, nkr);
  propagate_x(obj, idx);

  PROTECT(ret = allocMatrix(REALSXP, obj->nx, nd));
  qf_copy_x(obj, REAL(ret), nd, 0);
  UNPROTECT(1);
  return ret;
}


SEXP r_do_integrate(SEXP extPtr, SEXP vars, SEXP lambda, SEXP mu, 
		    SEXP drift, SEXP diffusion, SEXP nt, SEXP dt, 
		    SEXP padding) {
  quasse_fft *obj = (quasse_fft*)R_ExternalPtrAddr(extPtr);
  if ( obj == NULL )
    error("Corrupt QuaSSE integrator: ptr is NULL (are you using multicore?)");
  SEXP ret;
  int nkl = INTEGER(padding)[0], nkr = INTEGER(padding)[1];
  int ndat = LENGTH(lambda);
  double c_dt = REAL(dt)[0], c_nt = INTEGER(nt)[0];
  double *c_lambda=REAL(lambda), *c_mu=REAL(mu);
  double c_drift=REAL(drift)[0], c_diffusion=REAL(diffusion)[0];
  int i, idx, nd = LENGTH(vars) / obj->nx;
  
  idx = lookup(nd, obj->nd, obj->n_fft);
  if ( idx < 0 )
    error("Failed to find nd = %d\n", nd);

  qf_copy_x(obj, REAL(vars), nd, 1);

  obj->lambda = REAL(lambda);
  obj->mu = REAL(mu);
  for ( i = 0; i < ndat; i++ )
    obj->z[i] = exp(c_dt * (c_lambda[i] - c_mu[i]));

  qf_setup_kern(obj, c_drift, c_diffusion, c_dt, nkl, nkr);

  do_integrate(obj, c_nt, idx);

  obj->lambda = NULL;
  obj->mu = NULL;

  PROTECT(ret = allocMatrix(REALSXP, obj->nx, nd));
  qf_copy_x(obj, REAL(ret), nd, 0);
  UNPROTECT(1);
  
  return ret;
}

/*
  This one is much more specialised than r_do_integrate, and applies
  only to tips.  Assume that we get initial conditions in the order
    E, D[1], D[2], ..., D[nd_max]
  Now, we must do 'm' integrations with dimension
    nd_max, nd_max-1, nd_max-2, ..., nd_max-(m-1)
  where nd_max is the maximium possible extent given the integrator,
  and m is the number of different integrations, 1 <= m <= (nd_max-1)

  The D values must be ordered so that the lowest value D will be
  integrated for the longest, and the highest value D will be
  integrated for the shortest.
*/

SEXP r_do_tips(SEXP extPtr, SEXP vars, SEXP lambda, SEXP mu,
	       SEXP drift, SEXP diffusion, SEXP nt, SEXP dt, 
	       SEXP padding) {
  /* Setup directly copied from r_do_integrate, except that c_dt and
     c_nt are not initialised */
  quasse_fft *obj = (quasse_fft*)R_ExternalPtrAddr(extPtr);
  SEXP ret;
  int nkl = INTEGER(padding)[0], nkr = INTEGER(padding)[1];
  int ndat = LENGTH(lambda);
  double c_dt, c_nt;
  double *c_lambda=REAL(lambda), *c_mu=REAL(mu);
  double c_drift=REAL(drift)[0], c_diffusion=REAL(diffusion)[0];
  int i, idx;

  /* New setup */
  int nd, nx=obj->nx;
  int n_fft = obj->n_fft, n_fft_m1 = obj->n_fft - 1;

  if ( (LENGTH(vars) / obj->nx) != obj->nd[0] )
    error("Error 1\n");

  /* First; allocate space: All but the first cases will be nx * 2
     matrices, but the final one might be a matrix itself */
  PROTECT(ret = allocVector(VECSXP, n_fft));
  for ( i = 0; i < n_fft_m1; i++ )
    SET_VECTOR_ELT(ret, i, allocMatrix(REALSXP, nx, 2));
  SET_VECTOR_ELT(ret, n_fft_m1, 
		 allocMatrix(REALSXP, nx, obj->nd[n_fft_m1]));

  /* This bit proceeds exactly as r_do_integrate() */
  qf_copy_x(obj, REAL(vars), LENGTH(vars) / obj->nx, 1);

  obj->lambda = REAL(lambda);
  obj->mu = REAL(mu);

  /* New again */
  for ( idx = 0; idx < n_fft; idx++ ) {
    c_dt = REAL(dt)[idx];
    c_nt = INTEGER(nt)[idx];
    nd = obj->nd[idx];

    if ( c_nt > 0 ) {
      for ( i = 0; i < ndat; i++ )
	obj->z[i] = exp(c_dt * (c_lambda[i] - c_mu[i]));
      qf_setup_kern(obj, c_drift, c_diffusion, c_dt, nkl, nkr);
      do_integrate(obj, c_nt, idx);
    }

    if ( idx < (n_fft-1) )
      qf_copy_ED(obj, REAL(VECTOR_ELT(ret, idx)), nd-1);
    else
      qf_copy_x(obj, REAL(VECTOR_ELT(ret, idx)), nd, 0);
  }

  obj->lambda = NULL;
  obj->mu = NULL;
  
  UNPROTECT(1);
  
  return ret;
}

/* This does the memory allocation and plans the FFT transforms */
quasse_fft* make_quasse_fft(int n_fft, int nx, double dx, int *nd, 
			    int flags) {
  quasse_fft *obj = calloc(1, sizeof(quasse_fft));
  int ny = (((int)floor(nx/2)) + 1);
  int i, max_nd=1;
  for ( i = 0; i < n_fft; i++ )
    if ( nd[i] > max_nd )
      max_nd = nd[i];

  obj->n_fft = n_fft;
  obj->nx = nx;
  obj->ny = ny;
  obj->dx = dx;
  obj->nd = nd;

  obj->x   = fftw_malloc(max_nd *  nx    * sizeof(double));
  obj->y   = fftw_malloc(max_nd * (ny+1) * sizeof(complex));

  obj->z   = (double*)calloc(nx, sizeof(double));
  obj->wrk = (double*)calloc(nx, sizeof(double));

  obj->fft = (rfftw_plan_real**)calloc(n_fft, sizeof(rfftw_plan_real*));

  for ( i = 0; i < n_fft; i++ ) {
    obj->fft[i] = make_rfftw_plan_real(nd[i], nx, DIR_COLS,
				       obj->x, obj->y, flags);
  }
  
  /* Brownian kernel */
  obj->kern_x = fftw_malloc(nx     * sizeof(double));
  obj->kern_y = fftw_malloc((ny+1) * sizeof(complex));
  obj->kernel = make_rfftw_plan_real(1, nx, DIR_COLS, 
				     obj->kern_x, obj->kern_y, flags);
  
  return obj;
}

static void quasse_fft_finalize(SEXP extPtr) {
  quasse_fft *obj = (quasse_fft*)R_ExternalPtrAddr(extPtr);
  int i;
  /* Rprintf("Cleaning up\n"); */

  for ( i = 0; i < obj->n_fft; i++ ) {
    fftw_destroy_plan(obj->fft[i]->plan_f);
    fftw_destroy_plan(obj->fft[i]->plan_b);
  }
  free(obj->fft);
  free(obj->nd);

  fftw_free(obj->x);
  fftw_free(obj->y);

  free(obj->z);
  free(obj->wrk);

  fftw_destroy_plan(obj->kernel->plan_f);
  fftw_destroy_plan(obj->kernel->plan_b);

  fftw_free(obj->kern_x);
  fftw_free(obj->kern_y);

  free(obj);
}


void qf_copy_x(quasse_fft *obj, double *x, int nd, int copy_in) {
  int i, n = obj->nx * nd;
  double *fft_x = obj->x;
  if ( copy_in )
    for ( i = 0; i < n; i++ )
      fft_x[i] = x[i];
  else {
    for ( i = 0; i < n; i++ ) {
      x[i] = fft_x[i];
    }
  }
}

void qf_copy_ED(quasse_fft *obj, double *x, int idx) {
  int i, nx = obj->nx;
  double *fft_x = obj->x;
  for ( i = 0; i < nx; i++ )
    x[i] = fft_x[i];

  x += nx;
  fft_x = obj->x + idx*nx;

  for ( i = 0; i < nx; i++ )
    x[i] = fft_x[i];
}

void qf_setup_kern(quasse_fft *obj, double drift, double diffusion,
		   double dt, int nkl, int nkr) {
  const int nx = obj->nx;  
  int i;
  double x, *kern_x=obj->kern_x, tot=0, dx=obj->dx;
  double mean, sd;
  
  obj->nkl  = nkl;
  obj->nkr  = nkr;
  obj->npad = nkl + 1 + nkr;
  obj->ndat = nx - obj->npad;
  obj->drift = drift;
  obj->diffusion = diffusion;

  tot   = 0;
  mean  = - dt * drift;
  sd = sqrt(dt * diffusion);

  for ( i = 0, x = 0; i <= nkr; i++, x += dx )
    tot += kern_x[i] = dnorm(x, mean, sd, 0)*dx;
  for ( i = nkr + 1; i < nx - nkl; i++ )
    kern_x[i] = 0;
  for ( i = nx - nkl, x = -nkl*dx; i < nx; i++, x += dx )
    tot += kern_x[i] = dnorm(x, mean, sd, 0)*dx;

  for ( i = 0;        i <= nkr; i++ ) kern_x[i] /= tot;
  for ( i = nx - nkl; i < nx;   i++ ) kern_x[i] /= tot;

  fftw_execute(obj->kernel->plan_f);
}

void do_integrate(quasse_fft *obj, int nt, int idx) {
  int i, nkl=obj->nkl;
  for ( i = 0; i < nt; i++ ) {
    propagate_t(obj, idx);
    propagate_x(obj, idx);
    if ( ISNAN(obj->x[nkl]) )
      error("Integration failure at step %d\n", i);
  }
}

/* Lower level functions */
void propagate_t(quasse_fft *obj, int idx) {
  int ix, id, nx=obj->nx, ndat=obj->ndat, nd=obj->nd[idx];
  double *vars=obj->x, *d, *dd = obj->wrk;
  double e, tmp1, tmp2, lambda_x, mu_x, z_x;

  for ( ix = 0; ix < ndat; ix++ ) {
    lambda_x = obj->lambda[ix];
    mu_x     = obj->mu[ix];
    z_x      = obj->z[ix];
    e        = vars[ix];
    
    /* Update the E values */
    tmp1 = mu_x - lambda_x * e;
    tmp2 = z_x * (e - 1);
    vars[ix] = (tmp1 + tmp2 * mu_x) / (tmp1 + tmp2 * lambda_x);

    tmp1 = (lambda_x - mu_x) / 
      (z_x*lambda_x - mu_x + (1 - z_x)*lambda_x*e);
    /* Here is the D scaling factor */
    dd[ix] = z_x * tmp1 * tmp1;
  }

  /* Update the D values */
  for ( id = 1; id < nd; id++ ) {
    d = obj->x + nx * id;

    for ( ix = 0; ix < ndat; ix++ )
      if ( d[ix] < 0 )
	d[ix] = 0;
      else
	d[ix] *= dd[ix];
  }
}

void propagate_x(quasse_fft *obj, int idx) {
  double *x = obj->x, *wrk = obj->wrk;
  int i, id, nx = obj->nx;
  int nkl = obj->nkl, nkr = obj->nkr, npad = obj->npad;
  int nd=obj->nd[idx];

  /* TODO: I am not sure if these are flipped nkl/nkr */
  for ( i = 0; i < nkl; i++ )
    wrk[i] = x[i];
  for ( i = nx-npad-nkr; i < nx - npad; i++ )
    wrk[i] = x[i];

  convolve(obj->fft[idx], obj->kern_y);

  for ( i = 0; i < nkl; i++ )
    x[i] = wrk[i];
  for ( i = nx-npad-nkr; i < nx - npad; i++ )
    x[i] = wrk[i];

  /* Zeroing takes a little more work, now.  We might be able to get
     away with just zeroing the E though */
  for ( id = 0; id < nd; id++ ) {
    x = obj->x + (obj->nx)*(id+1) - npad;
    for ( i = 0; i < npad; i++ ) 
      x[i] = 0;
  }
}

void convolve(rfftw_plan_real *obj, complex *fy) {
  int nx = obj->nx, ny = obj->ny, nd = obj->nd, i, j;
  int nxd = nx * nd;
  double *x = obj->x;
  complex *y = obj->y;

  fftw_execute(obj->plan_f);

  for ( i = 0; i < nd; i++ )
    for ( j = 0; j < ny; j++, y++ )
      (*y) *= fy[j];

  fftw_execute(obj->plan_b);

  for ( i = 0; i < nxd; i++ )
    x[i] /= nx;
}

int lookup(int x, int *v, int len) {
  int i, idx=-1;
  for ( i = 0; i < len; i++ )
    if ( v[i] == x ) {
      idx = i;
      break;
    }

  return idx;
}

#endif
