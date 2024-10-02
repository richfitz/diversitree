/* Mkn - Markov k-state n-parameter character model */
#include <R.h>
#include "util.h"
#include "mkn-pij.h"

/* Simplified matrix multiplication, assuming straightforward sizes
   and zeroing the input.  GEMM does:
     Z = alpha X Y + beta Z
*/
void do_derivs_mknpij(int k, double *pars, const double *y, double *ydot) {
  double *Q = pars;
  do_gemm(Q, k, k, y, k, k, ydot);
}

void derivs_mknpij_gslode(int neqs, double t, double *pars, 
			  const double *y, double *dydt) {
  const int k = (int)sqrt(neqs);
  do_derivs_mknpij(k, pars, y, dydt);
}

/* TODO: Duplication within mkn-ode.c */
void initial_conditions_mknpij(int k, double *x_l, double *x_r, 
			       double *x_out) {
  int i;
  for ( i = 0; i < k; i++ )
    x_out[i] = x_l[i] * x_r[i];
}

/* 
   This assumes that everything is going to be stored in row major,
   rather than column major, order
*/
void mkn_core(int k, int n, int *order, int *children, double *pij,
	      double *branch_init, double *branch_base, double *lq) {
  int i, j, idx, idx_k;
  double *y_in, *y_out, q;

  for ( i = 0; i < n; i++ ) {
    idx = order[i];
    idx_k = idx * k;
    y_in = branch_init + idx_k;
    y_out = branch_base + idx_k;

    initial_conditions_mknpij(k, 
			      branch_base + k*children[idx*2],
			      branch_base + k*children[idx*2 + 1], 
			      y_in);

    do_gemm(pij + idx_k*k, k, k, y_in, k, 1, y_out);

    for ( q = 0.0, j = 0; j < k; j++ )
      q += y_out[j];

    for ( j = 0; j < k; j++ )
      y_out[j] /= q;

    lq[idx] = log(q);
  }

  /* Root */
  idx = order[n];
  idx_k = idx * k;
  y_in = branch_init + idx_k;
  initial_conditions_mknpij(k, 
			    branch_base + k*children[idx*2],
			    branch_base + k*children[idx*2 + 1], 
			    y_in);
}


void r_mkn_core(int *k, int *n, int *order, int *children, double *pij,
		double *branch_init, double *branch_base, double *lq) {
  mkn_core(*k, *n, order, children, pij, branch_init, branch_base, lq);
}

/* I probably need some of the cache:

   parent
   children
   root

   And the precomputed results:
   
   init
   base
   pij

   These need to be the *untransposed* calculations.
*/
SEXP r_asr_marginal_mkn(SEXP r_k, SEXP r_pars, SEXP r_nodes, 
			SEXP cache, SEXP res,
			SEXP root_f, SEXP rho) {
  const int n_states = INTEGER(r_k)[0];
  const int neq = n_states;
  int n_nodes = LENGTH(r_nodes), *nodes = INTEGER(r_nodes);

  /* I think these are the only elements of the cache that we need */
  int *parent   = INTEGER(VECTOR_ELT(cache, 0));
  int *children = INTEGER(VECTOR_ELT(cache, 1));
  int root      = INTEGER(VECTOR_ELT(cache, 2))[0];

  /* And these are the precomputed bits we need */
  double *r_init = REAL(VECTOR_ELT(res, 0));
  double *r_base = REAL(VECTOR_ELT(res, 1));
  double *r_lq   = REAL(VECTOR_ELT(res, 2));
  /* Spot 3 has 'vals' as of 0.9-2 */
  double *pij    = REAL(VECTOR_ELT(res, 4));
  int n_out = LENGTH(VECTOR_ELT(res, 2));

  /* These will be modified each time */
  double *lq   = (double*) R_alloc(n_out * neq, sizeof(double));
  double *init = (double*) R_alloc(n_out * neq, sizeof(double));
  double *base = (double*) R_alloc(n_out * neq, sizeof(double));
  /* And this is a pointer to the root variables within */
  double *root_vals = init + root * neq;

  SEXP ret, cpy_root_vals, cpy_lq, R_fcall, tmp;

  int idx, i, j, k;
  double *vals;

  if ( !Rf_isFunction(root_f) )
    Rf_error("root_f must be a function");
  if ( !Rf_isEnvironment(rho) )
    Rf_error("rho must be a function");

  PROTECT(ret = Rf_allocMatrix(REALSXP, n_states, n_nodes));
  PROTECT(cpy_root_vals = Rf_allocVector(REALSXP, neq));
  PROTECT(cpy_lq        = Rf_allocVector(REALSXP, n_out));

  for ( i = 0; i < n_nodes; i++ ) {
    idx = nodes[i];

    vals = REAL(ret) + n_states * i;

    for ( j = 0; j < n_states; j++ ) {
      /* Copy clean data back in */
      memcpy(lq,   r_lq,   n_out *       sizeof(double));
      memcpy(init, r_init, n_out * neq * sizeof(double));
      memcpy(base, r_base, n_out * neq * sizeof(double));

      for ( k = 0; k < n_states; k++ )
	if ( k != j )
	  init[neq * idx + k] = 0.0;

      asr_marginal_mkn_1(k, idx, root, parent, children, pij,
			 init, base, lq);

      memcpy(REAL(cpy_root_vals), root_vals, neq   * sizeof(double));
      memcpy(REAL(cpy_lq),        lq,        n_out * sizeof(double));
      PROTECT(R_fcall = Rf_lang4(root_f, r_pars, cpy_root_vals, cpy_lq));
      PROTECT(tmp = Rf_eval(R_fcall, rho));
      vals[j] = REAL(tmp)[0];
      UNPROTECT(2);
    }

    asr_normalise(n_states, vals);
  }

  UNPROTECT(3);
  return ret;
}

void asr_marginal_mkn_1(int k, int node, int root,
			int *parent, int *children,
			double *pij, double *branch_init, 
			double *branch_base, double *lq) {
  const int neq = k;
  int j, idx, idx_k, *kids;
  double *y_in, *y_out;
  double q;

  idx = node;
  idx_k = idx * k;
  y_in = branch_init + idx_k;
  y_out = branch_base + idx_k;

  while ( idx != root ) {
    /* Branch calculation */
    do_gemm(pij + idx_k*k, k, k, y_in, k, 1, y_out);

    /* Normalisation */
    q = 0.0;
    for ( j = 0; j < k; j++ )
      q += y_out[j];
    for ( j = 0; j < k; j++ )
      y_out[j] /= q;
    lq[idx] = log(q);

    /* Preparation for next iteration */
    idx = parent[idx];
    idx_k = idx * k;

    y_in = branch_init + idx_k;
    y_out = branch_base + idx_k;
    kids = children + idx * 2;
    initial_conditions_mknpij(neq, 
			      branch_base + neq*kids[0], 
			      branch_base + neq*kids[1], y_in);
  }
}

void asr_normalise(int n_states, double *vals) {
  int i;
  double maxp = R_NegInf, tot = 0.0;
  for ( i = 0; i < n_states; i++ )
    if (vals[i] > maxp)
      maxp = vals[i];
  for ( i = 0; i < n_states; i++ ) {
    if (ISNAN(vals[i]))
      vals[i] = 0.0;
    else
      tot += vals[i] = exp(vals[i] - maxp);
  }
  for ( i = 0; i < n_states; i++ )
    vals[i] /= tot;
}
