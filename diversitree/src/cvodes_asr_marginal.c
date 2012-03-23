#include <R.h>
#include <Rinternals.h>

#include "config.h"
#ifdef WITH_CVODES

#include "cvodes/include/cvodes/cvodes.h"
#include "cvodes/include/nvector/nvector_serial.h"  
#include "cvodes/include/sundials/sundials_types.h"

#include "cvodes_obj.h"
#include "cvodes_all_branches.h"

#include "mkn.h"

void asr_marginal_1(dt_obj *obj, double *pars, int node, 
		    int *parent);

/* Needed:
   extPtr: pointer containing the dt_obj object
   pars: Vector of parameters
   r_nodes: Vector of nodes at which we want ASR computed
   r_parent: Vector with parent relationships (see below)
   r_states: Vector indicating which of the equations represent states
   (logical or integer?)
   root_f: R function for computing the root.
   rho: Environment

   The 'anc' and 'parent' components are really important.
   The 'anc' component contains the node itself, followed by 
*/
SEXP r_asr_marginal(SEXP extPtr, SEXP r_pars, SEXP r_nodes,
		    SEXP r_states, SEXP r_parent,
		    SEXP root_f, SEXP rho) {
  dt_obj *obj = (dt_obj*)R_ExternalPtrAddr(extPtr);
  double *pars = REAL(r_pars);
  int neq = obj->neq, n_out = obj->n_out;
  int n_nodes = LENGTH(r_nodes), *nodes = INTEGER(r_nodes);
  int n_states = LENGTH(r_states), *states = INTEGER(r_states);
  int *parent = INTEGER(r_parent);

  int i, j, k, idx;
  double *lq, *init, *base, *vals, *root_vals;

  SEXP ret, r_root_vals, r_lq, R_fcall, tmp;

  lq   = (double*) R_alloc(n_out,       sizeof(double));
  init = (double*) R_alloc(n_out * neq, sizeof(double));
  base = (double*) R_alloc(n_out * neq, sizeof(double));

  memcpy(lq,   obj->lq,   n_out *       sizeof(double));
  memcpy(init, obj->init, n_out * neq * sizeof(double));
  memcpy(base, obj->base, n_out * neq * sizeof(double));

  /* The only argument checking */
  if ( !isFunction(root_f) )
    error("root_f must be a function");
  if ( !isEnvironment(rho) )
    error("rho must be a function");
  if ( LENGTH(r_pars) != obj->np )
    error("Incorrect length parameters.  Expected %d, got %d",
	  obj->np, LENGTH(r_pars));

  PROTECT(ret = allocMatrix(REALSXP, n_states, n_nodes));
  PROTECT(r_root_vals = allocVector(REALSXP, neq));
  PROTECT(r_lq        = allocVector(REALSXP, n_out));

  root_vals = obj->init + obj->root * neq;

  for ( i = 0; i < n_nodes; i++ ) {
    idx = nodes[i];
    
    vals = REAL(ret) + n_states * i;
    for ( j = 0; j < n_states; j++ ) {
      /* Copy clean data back in */
      memcpy(obj->lq,   lq,   n_out *       sizeof(double));
      memcpy(obj->init, init, n_out * neq * sizeof(double));
      memcpy(obj->base, base, n_out * neq * sizeof(double));

      for ( k = 0; k < n_states; k++ )
	if ( k != j )
	  obj->init[neq * idx + states[k]] = 0.0;

      asr_marginal_1(obj, pars, idx, parent);

      memcpy(REAL(r_root_vals), root_vals, neq   * sizeof(double));
      memcpy(REAL(r_lq),        obj->lq,   n_out * sizeof(double));
      PROTECT(R_fcall = lang4(root_f, r_pars, r_root_vals, r_lq));
      PROTECT(tmp = eval(R_fcall, rho));
      vals[j] = REAL(tmp)[0];
      UNPROTECT(2);
    }

    asr_normalise(n_states, vals);
  }
  
  UNPROTECT(3);
  return ret;
}

void asr_marginal_1(dt_obj *obj, double *pars, int node, 
		    int *parent) {
  RCvodesObj *integrator = obj->integrator;
  DtIcFun    ic = obj->ic;
  int root = obj->root, idx = node;
  double *init = obj->init, *base = obj->base, *lq = obj->lq;
  double *len = obj->len, *depth = obj->depth, eps = obj->eps;
  int *children = obj->children, *kids;
  int neq = obj->neq;
  double *vals = init + neq * idx;

  while ( idx != root ) {
    dt_cvodes_run(integrator, vals,
		  len + idx, 1, depth[idx], /* len[] and length(len) */
		  obj->comp_idx, obj->comp_n, eps,
		  &idx, base, lq);
    idx = parent[idx];
    vals = init + neq * idx;
    kids = children + idx * 2;
    ic(neq, base + neq*kids[0], base + neq*kids[1],
       pars, depth[idx], vals);
  }
}

#endif
