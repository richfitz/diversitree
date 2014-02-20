#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

/* #define VERBOSE */
#define SMKN_GROWTH_RATE 4
#define SMKN_MAX_SIZE    1000000
#define SCM_MAX_ATTEMPTS 100000
#include "scm-mkn.h"

static void smkn_info_finalize(SEXP extPtr);

SEXP r_smkn_alloc(SEXP k, SEXP n_hist) {
  smkn_info *obj;
  SEXP extPtr;
  obj = smkn_alloc(INTEGER(k)[0], INTEGER(n_hist)[0]);
  extPtr = R_MakeExternalPtr(obj, R_NilValue, R_NilValue);
  R_RegisterCFinalizer(extPtr, smkn_info_finalize);
  return extPtr;
}

SEXP r_smkn_set_pars(SEXP extPtr, SEXP pars) {
  smkn_info *obj = (smkn_info*)R_ExternalPtrAddr(extPtr);
  smkn_set_pars(obj, REAL(pars));
  return R_NilValue;
}

SEXP r_smkn_scm_run(SEXP extPtr, SEXP len,
		    SEXP state_beg, SEXP state_end, SEXP r_as01) {
  smkn_info *obj = (smkn_info*)R_ExternalPtrAddr(extPtr);
  SEXP ret = R_NilValue;
  int i, n_hist, as01 = LOGICAL(r_as01)[0];
  double *hist_t, *out_t, *out_s;
  int *hist_s;

  GetRNGstate(); 
  smkn_scm_run(obj, REAL(len)[0], INTEGER(state_beg)[0],
	       INTEGER(state_end)[0]);
  PutRNGstate(); /* below here is deterministic */

  n_hist = obj->n_hist;
  hist_t = obj->hist_t;
  hist_s = obj->hist_to;

  PROTECT(ret = allocMatrix(REALSXP, n_hist+1, 2));
  out_t = REAL(ret);
  out_s = REAL(ret) + (n_hist + 1);

  out_t[0] = 0.0;
  out_s[0] = as01 ? INTEGER(state_beg)[0] : INTEGER(state_beg)[0]+1;

  for ( i = 0; i < n_hist; i++ ) {
    out_t[i+1] = hist_t[i];
    out_s[i+1] = as01 ? hist_s[i] : hist_s[i]+1;
  }

  UNPROTECT(1);

  return ret;
}

SEXP r_smkn_scm_run_all(SEXP extPtr, SEXP pars, SEXP r_len, 
			SEXP r_state_beg, SEXP r_state_end, 
			SEXP r_as01, SEXP r_slim) {
  smkn_info *obj = (smkn_info*)R_ExternalPtrAddr(extPtr);
  SEXP ret, hist;
  int i, j, n_hist, n_hist_out, n = LENGTH(r_len);
  double *len = REAL(r_len);
  int *state_beg = INTEGER(r_state_beg), *state_end = INTEGER(r_state_end);
  int as01 = LOGICAL(r_as01)[0], slim = LOGICAL(r_slim)[0];
  double *hist_t, *out_t, *out_s;
  int *hist_s;

  if ( LENGTH(r_state_beg) != n )
    error("state_beg incorrect length");
  if ( LENGTH(r_state_end) != n )
    error("state_end incorrect length");

  smkn_set_pars(obj, REAL(pars));

  PROTECT(ret = allocVector(VECSXP, n));

  GetRNGstate();

  for ( i = 0; i < n; i++ ) {
    smkn_scm_run(obj, len[i], state_beg[i], state_end[i]);

    n_hist = obj->n_hist;
    n_hist_out = slim ? n_hist : n_hist + 1;

    hist_t = obj->hist_t;
    hist_s = obj->hist_to;
    
    /* List assignment protects the result so explicit
       PROTECT/UNPROTECT not neededed */
    SET_VECTOR_ELT(ret, i, allocMatrix(REALSXP, n_hist_out, 2));
    hist = VECTOR_ELT(ret, i);
    out_t = REAL(hist);
    out_s = REAL(hist) + n_hist_out;

    if ( slim ) {
      memcpy(out_t, hist_t, n_hist * sizeof(double));
      if ( as01 ) 
	for ( j = 0; j < n_hist; j++ ) out_s[j] = hist_s[j];
      else
	for ( j = 0; j < n_hist; j++ ) out_s[j] = hist_s[j] + 1;
    } else {
      out_t[0] = 0.0;
      out_s[0] = as01 ? state_beg[i] : state_beg[i]+1;

      memcpy(out_t + 1, hist_t, n_hist * sizeof(double));
      if ( as01 ) 
	for ( j = 0; j < n_hist; j++ ) out_s[j+1] = hist_s[j];
      else
	for ( j = 0; j < n_hist; j++ ) out_s[j+1] = hist_s[j] + 1;
    }
  }

  PutRNGstate();

  if ( slim )
    ret = smkn_slim(ret);

  UNPROTECT(1);
  return ret;
}

/* This function needs to be used carefully; In particular, the object
   returned and the original object 'obj' must not both be returned
   back to R */
SEXP smkn_slim(SEXP obj) {
  int i, j = 0, n = LENGTH(obj), nkeep = 0, *idx;
  SEXP r_idx, hist, ret;

  for ( i = 0; i < n; i++ )
    if ( nrows(VECTOR_ELT(obj, i)) > 0 )
      nkeep++;

  PROTECT(ret = allocVector(VECSXP, 2));
  PROTECT(r_idx = allocVector(INTSXP, nkeep));
  PROTECT(hist = allocVector(VECSXP, nkeep));

  idx = INTEGER(r_idx);
  for ( i = 0; i < n; i++ )
    if ( nrows(VECTOR_ELT(obj, i)) > 0 ) {
      idx[j] = i + 1;
      SET_VECTOR_ELT(hist, j, VECTOR_ELT(obj, i));
      j++;
    }

  SET_VECTOR_ELT(ret, 0, r_idx);
  SET_VECTOR_ELT(ret, 1, hist);
  UNPROTECT(3);

  return ret;
}

/* Functions required by _alloc() */
smkn_info* smkn_alloc(int k, int n_hist) {
  smkn_info *obj = Calloc(1, smkn_info);
  int np = k*(k - 1);
#ifdef VERBOSE
  Rprintf("Allocating with k = %d, n_hist = %d\n", k, n_hist);
#endif
  obj->k = k;

  obj->pars = Calloc(np,  double);
  obj->r    = Calloc(k,   double);
  obj->cp   = Calloc(np,  double);
  obj->perm = Calloc(np,  int);

  obj->n_hist = 0;
  obj->n_hist_max = n_hist;

  obj->hist_from = Calloc(n_hist, int);
  obj->hist_to   = Calloc(n_hist, int);

  obj->hist_t = Calloc(n_hist, double);

  return obj;
}

void smkn_cleanup(smkn_info *obj) {
#ifdef VERBOSE
  Rprintf("Cleaning permanantly sized objects\n");
#endif
  Free(obj->pars);
  Free(obj->r);
  Free(obj->cp);
  Free(obj->perm);

#ifdef VERBOSE
  Rprintf("Cleaning variably sized objects\n");
#endif
  Free(obj->hist_from);
  Free(obj->hist_to);
  Free(obj->hist_t);

  Free(obj);
}

static void smkn_info_finalize(SEXP extPtr) {
  smkn_cleanup((smkn_info*)R_ExternalPtrAddr(extPtr));
}

/* Function required by _set_pars() */
void smkn_set_pars(smkn_info* obj, double *pars) {
  int i, j, k = obj->k, k1 = obj->k - 1;
  int *perm, np = k * (k - 1);
  double *p, *cp, tot;

  memcpy(obj->pars, pars, np * sizeof(double));

  for ( i = 0; i < k; i++ ) {
    p    = obj->pars + i * k1;
    perm = obj->perm + i * k1;
    cp   = obj->cp   + i * k1;

    for ( j = 0; j < k1; j++ )
      cp[j] = p[j];

    for ( j = 0; j < k1; j++ )
      perm[j] = j;

    revsort(cp, perm, k1);
    
    for ( j = 1; j < k1; j++ )
      cp[j] += cp[j - 1];
    obj->r[i] = tot = cp[k1-1];
    for ( j = 0; j < k1; j++ )
      cp[j] /= tot;
  }
}

/* Functions required by _scm_run() */
int smkn_scm_run(smkn_info *obj, double len, 
		 int state_beg, int state_end) {
  int niter;

  for ( niter = 0; niter < SCM_MAX_ATTEMPTS; niter++ ) {
    smkn_sim(obj, state_beg, len);
    if ( obj->state == state_end )
      break;
  }

  if ( niter == SCM_MAX_ATTEMPTS )
    error("Realisation failed (too many attempts)");
  
  return 1;
}

double smkn_sim(smkn_info *obj, int x0, double len) {
  /* Temporary variables */
  double t = 0.0, dt;
  int state = x0, state_to;
#ifdef VERBOSE
  int niter = 0;
#endif

  smkn_init(obj, x0);

  for ( ;; ) {
    dt = rexp(1 / obj->r_tot);
    t = t + dt;

    if ( t > len ) {
      t = len;
      break;
    }

    state_to = smkn_pick_state(obj, state);
#ifdef VERBOSE
    Rprintf("%d (%2.5f): dt = %2.5f; %d -> %d \n", ++niter, t, dt,
	    state, state_to);
#endif

    smkn_evolve(obj, state, t, state_to);
    state = state_to;
  }

  return t;
}

void smkn_init(smkn_info *obj, int x0) {
  obj->state  = x0;
  obj->r_tot = obj->r[x0];
  obj->n_hist = 0;
}

int smkn_pick_state(smkn_info *obj, int state) {
  const int k = obj->k;
  int i, offset = (k-1)*state;
  double rU;
  double *cp = obj->cp + offset;
  int state_to;

  if ( k == 2 ) {
    state_to = state ? 0 : 1;
  } else {
    rU = unif_rand();
    for ( i = 0; i < k; i++ )
      if ( rU < cp[i] )
	break;

    state_to = obj->perm[offset + i];
    if ( state_to >= state )
      state_to++;
  }

  return state_to;
}

void smkn_evolve(smkn_info* obj, int state, double t, int state_to) {
  int i = obj->n_hist;  
  obj->state = state_to;
  obj->r_tot = obj->r[state_to];

  if ( obj->n_hist + 2 > obj->n_hist_max )
    smkn_grow_hist(obj);

  obj->hist_from[i] = state;
  obj->hist_to[i]   = state_to;
  obj->hist_t[i]    = t;
  obj->n_hist++;
}

void smkn_grow_hist(smkn_info *obj) {
  int n = obj->n_hist_max * SMKN_GROWTH_RATE;
#ifdef VERBOSE
  Rprintf("Growing hist from size %d to %d\n",
	  obj->n_hist_max, n);
#endif
  if ( n > SMKN_MAX_SIZE )
    error("Exceeding maximum allowed history size");

  obj->n_hist_max = n;

  obj->hist_from = Realloc(obj->hist_from, n, int);
  obj->hist_to   = Realloc(obj->hist_to,   n, int);
  obj->hist_t    = Realloc(obj->hist_t,    n, double);
}
