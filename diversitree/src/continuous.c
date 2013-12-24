#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "continuous.h"
#include "util.h"


static void dt_obj_cont_finalize(SEXP extPtr);
SEXP r_make_dt_obj_cont(SEXP cache, SEXP r_ic, SEXP r_br) {
  SEXP info = getListElement(cache, "info");
  int 
    neq = INTEGER(getListElement(info, "ny"))[0], 
    np  = INTEGER(getListElement(info, "np"))[0];
  /* Initial conditions and branches functions */
  DtIcFun ic = (DtIcFun) R_ExternalPtrAddr(r_ic);
  DtBrFun br = (DtBrFun) R_ExternalPtrAddr(r_br);

  /* Return object */
  dt_obj_cont *obj;
  SEXP extPtr;

  obj = (dt_obj_cont *)Calloc(1, dt_obj_cont);
  obj->neq  = neq;
  obj->n_out = LENGTH(getListElement(cache, "len"));
  obj->np   = np;
  obj->root = INTEGER(getListElement(cache, "root"))[0];

  obj->ic = ic;
  obj->br = br;

  /* Set up storage */
  obj->init = (double *)Calloc(obj->n_out * neq, double);
  obj->base = (double *)Calloc(obj->n_out * neq, double);
  obj->lq   = (double *)Calloc(obj->n_out,       double);

  /* Set up tips and internal branches */
  dt_cont_setup_tips(obj, cache);
  dt_cont_setup_internal(obj, cache);

  extPtr = R_MakeExternalPtr(obj, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(extPtr, dt_obj_cont_finalize);
  return extPtr;
}

/* Clean up

   NOTE: This is in approximately reverse order to the allocation
*/
static void dt_obj_cont_finalize(SEXP extPtr) {
  dt_obj_cont *obj = (dt_obj_cont*)R_ExternalPtrAddr(extPtr);

  Free(obj->init);
  Free(obj->base);
  Free(obj->lq);

  /* tips */
  Free(obj->tip_y);
  Free(obj->tip_len);
  Free(obj->tip_target);

  /* internals */
  Free(obj->order);
  Free(obj->children);
  Free(obj->len);
  Free(obj->depth);

  Free(obj);
}

void dt_cont_setup_tips(dt_obj_cont *obj, SEXP cache) {
  SEXP y, tip_y, tip_target_r;
  int i, idx, neq = obj->neq, *tip_target;
  int n_tip;

  y = getListElement(cache, "y");
  tip_y        = getListElement(y, "y");
  tip_target_r = getListElement(y, "target");
  tip_target   = INTEGER(tip_target_r);
  
  n_tip = obj->n_tip = LENGTH(tip_target_r);
  if ( nrows(tip_y) != neq || ncols(tip_y) != n_tip )
    error("Incorrect tip state dimensions");

  obj->tip_target = (int *)Calloc(n_tip, int);
  memcpy(obj->tip_target, tip_target, n_tip*sizeof(int));

  for ( i = 0; i < n_tip; i++ ) {
    idx = tip_target[i];
    memcpy(obj->init   + neq * idx,
	   REAL(tip_y) + neq * i,
	   neq * sizeof(double));
  }
}

void dt_cont_setup_internal(dt_obj_cont *obj, SEXP cache) {
  SEXP len, order, depth, children;
  int n_int, n_out;

  len      = getListElement(cache, "len");
  order    = getListElement(cache, "order");
  depth    = getListElement(cache, "depth");
  children = getListElement(cache, "children");

  n_out = obj->n_out;
  n_int = obj->n_int = LENGTH(order) - 1;

  obj->order    = (int    *)Calloc(n_int,   int);
  obj->children = (int    *)Calloc(n_out*2, int);
  obj->len      = (double *)Calloc(n_out,   double);
  obj->depth    = (double *)Calloc(n_out,   double);

  memcpy(obj->order,    INTEGER(order),      n_int*sizeof(int));
  memcpy(obj->children, INTEGER(children), 2*n_out*sizeof(int));
  memcpy(obj->len,      REAL(len),           n_out*sizeof(double));
  memcpy(obj->depth,    REAL(depth),         n_out*sizeof(double));
}

SEXP r_all_branches_cont(SEXP extPtr, SEXP r_pars) {
  dt_obj_cont *obj = (dt_obj_cont*)R_ExternalPtrAddr(extPtr);
  double *pars = REAL(r_pars);

  int neq, n_tip, n_int, *children, *order, *tip_target;
  double *len, *depth, *init, *base, *lq;

  DtIcFun    ic;
  DtBrFun    br;

  int i, idx, *kids;
  double tot = 0.0;
  SEXP ret, ret_vals;

  if ( obj == NULL )
    error("Corrupt pointer (are you using multicore?)");

  ic = obj->ic;
  br = obj->br;

  neq = obj->neq;
  n_tip = obj->n_tip;
  n_int = obj->n_int;
  tip_target = obj->tip_target;
  children = obj->children;
  order = obj->order;
  len = obj->len;
  depth = obj->depth;
  init = obj->init;
  base = obj->base;
  lq = obj->lq;

  if ( LENGTH(r_pars) != obj->np )
    error("Incorrect length parameters.  Expected %d, got %d",
	  obj->np, LENGTH(r_pars));

  for ( i = 0; i < n_tip; i++ ) {
    idx = tip_target[i];
    lq[idx] = br(init + neq * idx, len[idx], pars,
		 depth[idx], idx, base + neq * idx);
  }

  for ( i = 0; i < n_int; i++ ) {
    idx = order[i];
    kids = children + idx * 2;
    ic(neq, base + neq*kids[0], base + neq*kids[1],
       pars, depth[idx], init + neq * idx);

    lq[idx] = br(init + neq * idx, len[idx], pars,
		 depth[idx], idx, base + neq * idx);
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

SEXP r_get_vals_cont(SEXP extPtr) {
  dt_obj_cont *obj = (dt_obj_cont*)R_ExternalPtrAddr(extPtr);
  SEXP ret, r_base, r_init, lq;
  double *base, *init;
  int n_out = obj->n_out, neq = obj->neq;
  int i, idx;

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

  /* Strictly, we should set the root to NA */
  idx = obj->root * neq;
  for ( i = 0; i < neq; i++ )
    base[idx + i] = NA_REAL;

  UNPROTECT(4);
  return ret;
}

double branches_bm(double *vars_in, double len, double *pars, 
		   double t0, int idx, double *vars_out) {
  vars_out[0] = vars_in[0];
  vars_out[1] = vars_in[1] + pars[0] * len;
  vars_out[2] = 0.0;
  return vars_in[2];
}

/* Note use of expm1 here; saves numerical disaster when alpha is very
   very small (e.g. 1e-20); this is needed so that we converge
   gracefully on the alpha=0 case. */
double branches_ou_opt(double *vars_in, double len, double *pars, 
		       double t0, int idx, double *vars_out) {
  const double m = vars_in[0], v = vars_in[1], z = vars_in[2],
    s2 = pars[0], alpha = pars[1], theta = pars[2];

  if (alpha > 0) {
    vars_out[0] = exp(len * alpha) * (m - theta) + theta;
    vars_out[1] = expm1(2*len*alpha) * s2 / (2*alpha) +
      exp(2*len*alpha) * v;
  } else {
    vars_out[0] = m;
    vars_out[1] = v + len * s2;
  }
  vars_out[2] = 0.0;

  return len * alpha + z;
}

double branches_ou_noopt(double *vars_in, double len, double *pars, 
		   double t0, int idx, double *vars_out) {
  const double m = vars_in[0], v = vars_in[1], z = vars_in[2],
    s2 = pars[0], alpha = pars[1];

  if (alpha > 0) {
    vars_out[0] = m;
    vars_out[1] = v - s2 *
      exp(-2 * alpha * t0) * expm1(-2 * alpha * len) / (2 * alpha);
  } else {
    vars_out[0] = m;
    vars_out[1] = v + len * s2;
  }
  vars_out[2] = 0.0;

  return z;
}

double branches_eb(double *vars_in, double len, double *pars,
		  double t0, int idx, double *vars_out) {
  const double m = vars_in[0], v = vars_in[1], z = vars_in[2],
    sigma2 = pars[0], a = pars[1], tr = pars[2];
  double s0, s1;
  if (a != 0) {
    s1 = tr - t0;
    s0 = s1 - len;
    len = (exp(a * s1) - exp(a * s0))/a;
  }
  vars_out[0] = m;
  vars_out[1] = v + sigma2 * len;
  vars_out[2] = 0.0;
  return z;
}

double branches_lambda(double *vars_in, double len, double *pars,
		       double t0, int idx, double *vars_out) {
  const double m = vars_in[0], v = vars_in[1], z = vars_in[2];
  const double sigma2 = pars[0], lambda = pars[1], tr = pars[2];
  const int n_tip = (int)pars[3];
  double len_scaled = len * lambda;
  if (idx <= n_tip)
    len_scaled += (1 - lambda) * (tr - t0);

  vars_out[0] = m;
  vars_out[1] = v + sigma2 * len_scaled;
  vars_out[2] = 0.0;

  return z;
}

/* Shared between bm and ou */
void initial_conditions_bm(int neq, double *vars_l, double *vars_r,
			   double *pars, double t, 
			   double *vars_out) {
  double m_l, m_r, v_l, v_r, vv;
  m_l = vars_l[0];
  m_r = vars_r[0];
  v_l = vars_l[1];
  v_r = vars_r[1];
  vv  = v_l + v_r;

  vars_out[0] = (m_l * v_r + m_r * v_l) / vv;
  vars_out[1] = (v_l * v_r) / vv;
  vars_out[2] = -(m_l - m_r)*(m_l - m_r) / (2 * vv) - 
    log(2 * M_PI * vv) / 2;
}

/* PGLS: tip_y is cache$y$y, but in tip order because of the
   reordering done by this point already).  Also, this code duplicates
   code from dt_cont_setup_tips.

   TODO: Wrap this up nicely (look in model-pgls.R for use).
 */
SEXP r_dt_cont_reset_tips(SEXP extPtr, SEXP tip_y) {
  dt_obj_cont *obj = (dt_obj_cont*)R_ExternalPtrAddr(extPtr);
  double *y = REAL(tip_y);
  int i, idx, neq = obj->neq, n_tip = obj->n_tip;
  if (LENGTH(tip_y) != neq * n_tip)
    error("Wrong length tip_y - expected %d, got %d",
	  neq * n_tip, LENGTH(tip_y));

  for (i = 0; i < n_tip; i++) {
    idx = obj->tip_target[i];
    memcpy(obj->init   + neq * idx,
	   y           + neq * i,
	   neq * sizeof(double));
  }

  return R_NilValue;
}
