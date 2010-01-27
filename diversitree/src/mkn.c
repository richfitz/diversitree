/* Mkn - Markov k-state n-parameter character model */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Rdynload.h>

static double *parms_mkn;
void initmod_mkn(void (* odeparms)(int *, double *)) {
  /* TODO: I should check here about the lengths of parameters, but I
     won't bother; it is not clear how best to do this, anyway.  Most
     of the checking should be done in the R end of things.  Because
     this is going to get an R object, it should be much easier to the
     checking there. */
  DL_FUNC get_deSolve_gparms = 
    R_GetCCallable("deSolve", "get_deSolve_gparms");
  parms_mkn = REAL(get_deSolve_gparms());
} 

/* Simplified matrix multiplication, assuming straightforward sizes
   and zeroing the input.  GEMM does:
     Z = alpha X Y + beta Z
*/
void do_gemm(double *x, int nrx, int ncx,
             double *y, int nry, int ncy,
             double *z) {
  const char *trans = "N";
  double one = 1.0, zero = 0.0;
  F77_CALL(dgemm)(trans, trans, &nrx, &ncy, &ncx, &one,
                  x, &nrx, y, &nry, &zero, z, &nrx);
}

void derivs_mkn(int *neq, double *t, double *y, double *ydot,
		double *yout, int *ip) {
  const int k = *neq;
  do_gemm(parms_mkn, k, k, y, k, 1, ydot);
}

void derivs_mkn_pij(int *neq, double *t, double *y, double *ydot,
		    double *yout, int *ip) {
  const int k = (int)sqrt(*neq);
  do_gemm(parms_mkn, k, k, y, k, k, ydot);
}

void r_gemm(double *x, int *nrx, int *ncx,
            double *y, int *nry, int *ncy,
            double *z) {
  do_gemm(x, *nrx, *ncx, y, *nry, *ncy, z);
}

void initial_conditions_mkn(int k, double *x_l, double *x_r, 
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

    initial_conditions_mkn(k, 
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
  initial_conditions_mkn(k, 
			 branch_base + k*children[idx*2],
			 branch_base + k*children[idx*2 + 1], 
			 y_in);
}


void r_mkn_core(int *k, int *n, int *order, int *children, double *pij,
		double *branch_init, double *branch_base, double *lq) {
  mkn_core(*k, *n, order, children, pij, branch_init, branch_base, lq);
}
