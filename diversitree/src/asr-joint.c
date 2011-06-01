#include <R.h>
#include <Rinternals.h>

/* Ripped from the R code; I can probably use this to speed things up
   later? - If I drop the permutations or save them in the main
   section this will be a bit faster... */
int ProbSampleOne_tmp(int n, double *p, int *perm) {
  double rU, tot = 0;
  int i, j;
  int nm1 = n - 1;

  /* record element identities */
  for (i = 0; i < n; i++)
    perm[i] = i;

  /* sort the probabilities into descending order */
  revsort(p, perm, n);

  /* compute cumulative probabilities */
  for (i = 1 ; i < n; i++)
    p[i] += p[i - 1];
  tot = p[i - 1];
  for (i = 0; i < n; i++)
    p[i] /= tot;

  /* compute the sample */
  rU = unif_rand();
  for (j = 0; j < nm1; j++)
    if (rU <= p[j])
      break;

  return perm[j];
}

SEXP r_sample(SEXP r_root_p) {
  int k = LENGTH(r_root_p);
  int *perm = (int*)R_alloc(k, sizeof(int));
  double *pr = (double*)R_alloc(k, sizeof(double));
  int ret;

  GetRNGstate();
  memcpy(pr, REAL(r_root_p), k * sizeof(double));
  ret = ProbSampleOne_tmp(k, pr, perm);
  PutRNGstate();
  
  return ScalarInteger(ret);
}


/* 
   li  must be k   x (n_tip * (n_tip - 1)) (transpose of the default)
   pij must be k*k x (n_tip * (n_tip - 1))

 */
SEXP r_do_asr_joint(SEXP r_k, SEXP r_order, SEXP r_parent,
		    SEXP r_li, SEXP r_pij, SEXP r_root_p, SEXP r_as_01) {
  const int k = INTEGER(r_k)[0], len = LENGTH(r_order);
  const int n_tip = len + 1, as_01 = LOGICAL(r_as_01)[0];
  const int k2 = k * k;
  int *order = INTEGER(r_order), *parent = INTEGER(r_parent);
  SEXP ret;
  int *states;
  int *perm = (int*)R_alloc(k, sizeof(int));
  double *pr = (double*)R_alloc(k, sizeof(double));
  double *li = REAL(r_li), *pij = REAL(r_pij), *li_i, *pij_i;
  int idx, i, j, l;

  GetRNGstate();
  PROTECT(ret = allocVector(INTSXP, len));
  states = INTEGER(ret);

  /* Sample root */
  memcpy(pr, REAL(r_root_p), k * sizeof(double));
  states[0] = ProbSampleOne_tmp(k, pr, perm);

  for ( i = 1; i < len; i++ ) {
    idx = order[i];
    j = states[parent[idx] - n_tip];
    li_i  = li  + idx * k;
    pij_i = pij + idx * k2;
    for ( l = 0; l < k; l++ )
      pr[l] = li_i[l] * pij_i[l*k+j];
    states[idx - n_tip] = ProbSampleOne_tmp(k, pr, perm);
  }

  if ( !as_01 )
    for ( i = 0; i < len; i++ )
      states[i]++;

  PutRNGstate();
  UNPROTECT(1);

  return ret;
}
