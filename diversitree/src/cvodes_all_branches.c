#include <R.h>
#include <Rinternals.h>

#include "config.h"
#ifdef WITH_CVODES

#include "cvodes/include/cvodes/cvodes.h"
#include "cvodes/include/nvector/nvector_serial.h"  
#include "cvodes/include/sundials/sundials_types.h"

#include "cvodes_obj.h"
#include "cvodes_all_branches.h"

/* TODO: This does not deal with presets at all, which means that it
   cannot be used with the unresolved clade method, or with the
   "split" methods */

/* TODO: Get rid of this, or set to zero */
/* Basically, we will flake out and not run the integration of a
   branch length is too small */
#define SMALL_STEP 1e-14

/* Maximum recusion depth.  I doubt that we would ever get there, but
   if we fail at 10, we're going to fail anyway (this breaks the
   branch into up to 2^10 = 1024 segments */
#define MAX_DEPTH 10


/* As usual, functions preceeded by 'r_' are intended/able to be
   called from R */

/* Create a 'dt_obj' object that contains all the necessary things for
   running the tree calculations.

   cache: A C-ized cache
   neq:   Number of equations that we will be solving
   np:    Number of parameters
   rhs:   Address of the rhs (derivative) function
   ic:    Address of the initial conditions function
   rtol:  Scalar relative tolerance
   atol:  Vector (length neq) absolute tolerance 
*/
static void dt_obj_finalize(SEXP extPtr);
SEXP r_make_dt_obj(SEXP cache, SEXP r_neq, SEXP r_np,
		   SEXP r_rhs, SEXP r_ic, SEXP r_rtol, SEXP r_atol,
		   SEXP r_eps) {
  int neq = INTEGER(r_neq)[0], np = INTEGER(r_np)[0];
  CVRhsFn rhs = (CVRhsFn) R_ExternalPtrAddr(r_rhs);
  DtIcFun ic  = (DtIcFun) R_ExternalPtrAddr(r_ic);
  dt_obj *obj;
  SEXP extPtr, comp_idx;

  if ( LENGTH(r_atol) != neq )
    error("Incorrect length of atol");

  obj = calloc(1, sizeof(dt_obj));
  obj->neq  = neq;
  obj->np   = np;

  obj->root = INTEGER(getListElement(cache, "root"))[0];

  /* 1: build integrator */
  obj->integrator = make_cvodes(neq, np, rhs, 
				REAL(r_rtol)[0], REAL(r_atol));
  if ( obj->integrator == NULL )
    error("Failed to build integrator");

  /* 2: Set up initial conditions */
  obj->ic = ic;

  /* 3: Set up under flow compensation */
  comp_idx = getListElement(cache, "comp.idx");
  obj->comp_n   = LENGTH(comp_idx);
  obj->comp_idx = calloc(obj->comp_n, sizeof(int));
  memcpy(obj->comp_idx, INTEGER(comp_idx), obj->comp_n*sizeof(int));
  obj->eps = REAL(r_eps)[0];

  /* 4: Set up tips and internal branches */
  dt_setup_tips(obj, cache);
  dt_setup_internal(obj, cache);

  /* 5: Set up storage */
  obj->init = calloc(obj->n_out * neq, sizeof(double));
  obj->base = calloc(obj->n_out * neq, sizeof(double));
  obj->lq   = calloc(obj->n_out,       sizeof(double));

  extPtr = R_MakeExternalPtr(obj, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(extPtr, dt_obj_finalize);
  return extPtr;
}

/* Extract the values that would normally go for intermediate
   calculation.  This is awkward as we have to transpose the init and
   base arrays
*/
SEXP r_get_vals(SEXP extPtr) {
  dt_obj *obj = (dt_obj*)R_ExternalPtrAddr(extPtr);
  SEXP ret, r_base, r_init, lq;
  double *base, *init;
  int n_out = obj->n_out, neq = obj->neq;
  int i, j, k;
  int tip_n, *tip_target;
  double *tip_y;

  PROTECT(ret = allocVector(VECSXP, 3));
  PROTECT(r_init = allocMatrix(REALSXP, neq, n_out));
  PROTECT(r_base = allocMatrix(REALSXP, neq, n_out));

  PROTECT(lq   = allocVector(REALSXP, n_out));
  SET_VECTOR_ELT(ret, 0, r_init);
  SET_VECTOR_ELT(ret, 1, r_base);
  SET_VECTOR_ELT(ret, 2, lq);

  base = REAL(r_base);
  init = REAL(r_init);

  memcpy(base,     obj->base, n_out * neq * sizeof(double));
  memcpy(init,     obj->init, n_out * neq * sizeof(double));
  memcpy(REAL(lq), obj->lq,   n_out *       sizeof(double));

  /* And then we have to copy the tip information in too... */
  for ( i = 0; i < obj->tip_types; i++ ) {
    tip_n      = obj->tip_n[i];
    tip_y      = obj->tip_y[i];
    tip_target = obj->tip_target[i];
    for ( j = 0; j < neq; j++ )
      for ( k = 0; k < tip_n; k++ ) 
	memcpy(init + tip_target[k] * neq, tip_y, neq * sizeof(double));
  }

  /* Strictly, we should set the root to NA */
  j = obj->root * neq;
  for ( i = 0; i < neq; i++ )
    base[j + i] = NA_REAL;

  UNPROTECT(4);
  return ret;
}

/* Run the 'all_branches()' calculation, similarly to the R verision.
   As everything has been set up correctly, we need only the
   parameters */
SEXP r_all_branches(SEXP extPtr, SEXP r_pars) {
  dt_obj *obj = (dt_obj*)R_ExternalPtrAddr(extPtr);
  double *pars = REAL(r_pars);

  RCvodesObj *integrator;
  DtIcFun    ic;
  int neq, tip_types, n_int, *children, *order;
  double eps, *len, *depth, *init, *base, *lq;

  int i, idx, *kids;
  double t0 = 0.0, tot = 0.0;
  SEXP ret, ret_vals;

  /* TODO:
     Ideally what we would do here is check *before* using, and
     recreate the object if needed.  This will require an additional
     .Call() and if statement each time.
   */
  if ( obj == NULL )
    error("Corrupt pointer (are you using multicore?)");

  integrator = obj->integrator;
  ic         = obj->ic;

  neq = obj->neq;
  tip_types = obj->tip_types;
  n_int = obj->n_int;
  children = obj->children;
  order = obj->order;
  len = obj->len;
  depth = obj->depth;
  init = obj->init;
  base = obj->base;
  lq = obj->lq;
  eps = obj->eps;

  /* cvodes_set_pars copies the parameter vector; no chance of back
     writing */
  if ( LENGTH(r_pars) != obj->np )
    error("Incorrect length parameters.  Expected %d, got %d",
	  obj->np, LENGTH(r_pars));
  cvodes_set_pars(integrator, pars);

  for ( i = 0; i < tip_types; i++ ) {
    dt_cvodes_run(integrator, obj->tip_y[i], obj->tip_len[i],
		  obj->tip_n[i], t0,
		  obj->comp_idx, obj->comp_n, eps,
		  obj->tip_target[i], base, lq);
  }

  for ( i = 0; i < n_int; i++ ) {
    idx = order[i];
    kids = children + idx * 2;
    ic(neq, base + neq*kids[0], base + neq*kids[1],
       pars, depth[idx], init + neq * idx);

    dt_cvodes_run(integrator, init + neq * idx,
		  len + idx, 1, depth[idx],
		  obj->comp_idx, obj->comp_n, eps,
		  &idx, base, lq);
  }

  idx = obj->root;
  kids = children + idx * 2;
  ic(neq, base + neq*kids[0], base + neq*kids[1],
     pars, depth[idx], init + neq * idx);

  lq[obj->root] = 0.0; /* defensive */
  for ( i = 0; i < obj->n_out; i++ )
    tot += lq[i];

  PROTECT(ret = allocVector(VECSXP, 2));
  PROTECT(ret_vals = allocVector(REALSXP, neq));
  SET_VECTOR_ELT(ret, 0, ScalarReal(tot));
  SET_VECTOR_ELT(ret, 1, ret_vals);
  memcpy(REAL(ret_vals), init + obj->root * neq, neq*sizeof(double));
  UNPROTECT(2);
    
  return ret;
}

/* Support for setting up the dt_obj object */
void dt_setup_tips(dt_obj *obj, SEXP cache) {
  SEXP y, el, tip_y, tip_len, tip_target;
  int i, n_i, neq = obj->neq;

  y = getListElement(cache, "y");
  obj->tip_types  = LENGTH(y);
  obj->tip_n      = calloc(obj->tip_types, sizeof(int));
  obj->tip_y      = calloc(obj->tip_types, sizeof(double*));
  obj->tip_len    = calloc(obj->tip_types, sizeof(double*));
  obj->tip_target = calloc(obj->tip_types, sizeof(int*));

  for ( i = 0; i < obj->tip_types; i++ ) {
    el = VECTOR_ELT(y, i);

    tip_y      = getListElement(el, "tip.y");
    tip_len    = getListElement(el, "tip.len");
    tip_target = getListElement(el, "tip.target");
    n_i = LENGTH(tip_target);

    obj->tip_n[i]      = n_i;
    obj->tip_y[i]      = calloc(neq, sizeof(double));
    obj->tip_len[i]    = calloc(n_i, sizeof(double));
    obj->tip_target[i] = calloc(n_i, sizeof(int));

    if ( LENGTH(tip_y) != neq )
      error("Incorrect size initial conditions");

    memcpy(obj->tip_y[i],   REAL(tip_y),   neq*sizeof(double));
    memcpy(obj->tip_len[i], REAL(tip_len), n_i*sizeof(double));
    memcpy(obj->tip_target[i], INTEGER(tip_target), n_i*sizeof(int));
  }
}

void dt_setup_internal(dt_obj *obj, SEXP cache) {
  SEXP len, order, depth, children;
  int n_int, n_out;

  len      = getListElement(cache, "len");
  order    = getListElement(cache, "order");
  depth    = getListElement(cache, "depth");
  children = getListElement(cache, "children");

  n_out = obj->n_out = LENGTH(len);
  n_int = obj->n_int = LENGTH(order) - 1;

  obj->order    = calloc(n_int,   sizeof(int));
  obj->children = calloc(n_out*2, sizeof(int));
  obj->len      = calloc(n_out,   sizeof(double));
  obj->depth    = calloc(n_out,   sizeof(double));

  memcpy(obj->order,    INTEGER(order),      n_int*sizeof(int));
  memcpy(obj->children, INTEGER(children), 2*n_out*sizeof(int));
  memcpy(obj->len,      REAL(len),           n_out*sizeof(double));
  memcpy(obj->depth,    REAL(depth),         n_out*sizeof(double));
}

/* Clean up

   NOTE: This is in approximately reverse order to the allocation
*/
static void dt_obj_finalize(SEXP extPtr) {
  dt_obj *obj = (dt_obj*)R_ExternalPtrAddr(extPtr);
  int i;

  free(obj->init);
  free(obj->base);
  free(obj->lq);

  for ( i = 0; i < obj->tip_types; i++ ) {
    free(obj->tip_y[i]);
    free(obj->tip_len[i]);
    free(obj->tip_target[i]);
  }
  free(obj->tip_n);
  free(obj->tip_y);
  free(obj->tip_len);
  free(obj->tip_target);

  free(obj->order);
  free(obj->children);
  free(obj->len);
  free(obj->depth);

  free(obj->comp_idx);

  cvodes_cleanup(obj->integrator);
  
  free(obj);
}

/* Utility functions that control the running */
void dt_cvodes_run(RCvodesObj *obj, double *y_in,
		   double *len, int nt, double t0,
		   int *comp_idx, int comp_n, double eps,
		   int *target, double *base, double *lq) {
  int neq = obj->neq;
  double t = t0, t1, *ydat = NV_DATA_S(obj->y), *out, *y;
  int i, j;
  double tmp, lq_cum = 0.0;
  memcpy(ydat, y_in, neq * sizeof(double));

  for ( i = 0; i < nt; i++ ) {
    j = target[i];
    out = base + neq*j;
    t1 = t0 + len[i];

    if ( t1 - t > SMALL_STEP ) {
      tmp = dt_cvodes_run_1(obj, t, t1, comp_idx, comp_n, eps);
      if ( tmp == R_NegInf ) {
	y = (i == 0) ? y_in : base + neq * target[i-1];
	tmp = dt_cvodes_run_multi(obj, y, t, t1, comp_idx, comp_n, 
				  eps, 1);
      }
      lq_cum += tmp;
    }

    memcpy(out, ydat, neq * sizeof(double));
    lq[j] = lq_cum;

    t = t1;
  }
}

double dt_cvodes_run_1(RCvodesObj *obj, double t0, double t1, 
		       int *comp_idx, int comp_n, double eps) {
  int i, flag, ok=1;
  void *cvode_mem = obj->cvode_mem;
  double tmp, q, lq = 0, t, *ydat = NV_DATA_S(obj->y);

  CVodeReInit(obj->cvode_mem, t0, obj->y);
  flag = CVode(cvode_mem, t1, obj->y, &t, CV_NORMAL);
  if ( flag != CV_TOO_CLOSE && cvodes_check_flag(&flag, "CVode", 1) )
    error("Solver failed (message probably above)");

  if ( comp_n > 0 ) {
    q = 0;
    for ( i = 0; i < comp_n; i++ ) {
      tmp = ydat[comp_idx[i]];
      if ( tmp < 0 ) 
	ok = 0;
      q += tmp;
    }

    if ( ok && q > eps ) {
      lq = log(q);
      for ( i = 0; i < comp_n; i++ )
	ydat[comp_idx[i]] /= q;
    } else {
      /* This will flag a recoverable failure */
      return R_NegInf;
    }
  }

  return lq;
}

double dt_cvodes_run_multi(RCvodesObj *obj, double *y, 
			   double t0, double t1, 
			   int *comp_idx, int comp_n, double eps,
			   int depth) {
  double *ydat = NV_DATA_S(obj->y), t12 = (t0 + t1)/2, lq1, lq2;
  if ( depth > MAX_DEPTH )
    error("Recursion limit reached");

  memcpy(ydat, y, obj->neq * sizeof(double));

  lq1 = dt_cvodes_run_1(obj, t0,  t12, comp_idx, comp_n, eps);
  if ( lq1 == R_NegInf )
    lq1 = dt_cvodes_run_multi(obj, y, t0, t12, comp_idx, comp_n, eps,
			      depth+1);
  lq2 = dt_cvodes_run_1(obj, t12, t1,  comp_idx, comp_n, eps);
  if ( lq2 == R_NegInf )
    lq2 = dt_cvodes_run_multi(obj, y, t0, t12, comp_idx, comp_n, eps,
			      depth+1);

  return lq1 + lq2;
}

/* R interface to the above functions for testing... */
SEXP r_dt_cvodes_run(SEXP extPtr, SEXP y, SEXP len, SEXP pars, 
		     SEXP t0, SEXP idx, SEXP eps, SEXP target) {
  RCvodesObj *obj = (RCvodesObj*)R_ExternalPtrAddr(extPtr);
  int nt = LENGTH(len), neq = obj->neq;
  SEXP ret, ret_lq, ret_base;

  PROTECT(ret = allocVector(VECSXP, 2));
  PROTECT(ret_lq = allocVector(REALSXP, nt));
  PROTECT(ret_base = allocMatrix(REALSXP, neq, nt));
  SET_VECTOR_ELT(ret, 0, ret_lq);
  SET_VECTOR_ELT(ret, 1, ret_base);

  cvodes_set_pars(obj, REAL(pars));
  dt_cvodes_run(obj, REAL(y), 
		REAL(len), LENGTH(len), REAL(t0)[0],
		INTEGER(idx), LENGTH(idx), REAL(eps)[0],
		INTEGER(target), REAL(ret_base), REAL(ret_lq));

  UNPROTECT(3);
  return ret;
}

/* Utility function for accessing list elements by name.  This is
   needed to stop the argument list getting out of control */
SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i; 
  for ( i = 0; i < length(list); i++ )
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) { 
      elmt = VECTOR_ELT(list, i); 
      break; 
    }

  if ( elmt == R_NilValue )
    error("%s missing from list", str);

  return elmt;
} 

/* I wonder if the problem is the tips? */
SEXP r_dt_debug(SEXP extPtr) {
  dt_obj *obj = (dt_obj*)R_ExternalPtrAddr(extPtr);
  SEXP ret;
  int tip_types = obj->tip_types, neq = obj->neq;
  double *tmp;
  int i;

  PROTECT(ret = allocMatrix(REALSXP, neq, tip_types));
  for ( i = 0; i < tip_types; i++ ) {
    tmp = REAL(ret) + i*neq;
    memcpy(tmp, obj->tip_y[i], neq*sizeof(double));
  }

  UNPROTECT(1);

  return ret;
}

#endif
