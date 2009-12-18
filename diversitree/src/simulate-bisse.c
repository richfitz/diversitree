#include <R.h>
#include <Rmath.h>
/*
#include <Rinternals.h>
*/

/* Ripped from the R code; I can probably use this to speed things up
   later? - If I drop the permutations or save them in the main
   section this will be a bit faster... */
int ProbSampleOne(int n, double *p, int *perm) {
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

int SampleOne(int n) {
  return (int)(n * unif_rand());
}

/* This should work */
/*
void r_ProbSampleOne(int *n, double *p, int *ans) {
  int *perm = (int*)R_alloc(*n, sizeof(int));
  GetRNGstate();
  *ans = ProbSampleOne(*n, p, perm);
  PutRNGstate();
}
void r_SampleOne(int *n, int *ans) {
  GetRNGstate();
  *ans = SampleOne(*n);
  PutRNGstate();
}
*/

/*
  'state' is defined by the elements
  idx: a vector of indices
  parent: parent index
  state: character state
  extinct: is the vector extinct?
  split: has the vector split?
  len: length of a segment

  There are bunch of problems here with pausing/resuming I think.
  To get this done efficiently, I might need to work in a
  .Call/PROTECT setting?
 */
/*
  It is simple to put this into an infinite loop if starting in a
  state that has no speciation rate with max_t = Inf.  I have not
  checked for this anywhere.
 */
int simulate_bisse(double *pars, int max_taxa, double max_t, 
		   int *parent, int *states, int *extinct, 
		   int *split, double *start, double *len,
		   int hist[][3], double *hist_t,
		   int n_entries, int n, double *t_start) {
  int *lineages;
  int i, j, k, n_taxa, lineage, state, type, n_hist=0;
  int n_i[2], perm[3];
  double r_i[2], r_n[2], r_tot, t=*t_start, dt;
  double p[3];

  /* 
     First, process the parameters; these come in the order
       l0, m0, q01,  l1, m1, q10
     WARNING: This is different to the R version!
     r_i gives the per-lineage rates over all events.  If I store the
     cumulative probabilities over different events, this would be the
     place that I would do it.
  */
  r_i[0] = pars[0] + pars[1] + pars[2];
  r_i[1] = pars[3] + pars[4] + pars[5];

  /* 
     Next, determine which lineages are "current" in the given data;
     these are lineages that have not gone extinct or split.  At the
     same time, count how many species are in each state (n_i) and how
     many species we have total (n_taxa).
  */
  lineages = (int*)R_alloc(n, sizeof(int));
  n_i[0] = n_i[1] = 0;
  for ( i = 0, j = 0; i < n_entries; i++ ) {
    if ( !extinct[i] && !split[i] ) {
      lineages[j++] = i;
      (n_i[states[i]])++;
    }
  }
  n_taxa = n_i[0] + n_i[1];

  /* This is the main loop */
  while ( n_taxa <= max_taxa && n_taxa > 0 ) {
    /* Need to check that we can fit additional speciation events
     *before* doing any RNG, so that results will be the same as n
     changes.  I could possibly do this _after_ a speciation event and
     avoid a few comparisons, but this is probably safter to guarantee
     no bias (if we bailed only after shorter intervals or speciation
     events, this might bias trees towards fewer late events and fewer
     speciation events without some care) */
    if ( n_entries + 2 > n ) {
      /*error("History has got too big; I will deal with this later.");*/
      *t_start = -t;
      return n_entries;
    }
      
    /* Starting round */
    for ( i = 0; i < 2; i++ )
      r_n[i] = n_i[i] * r_i[i];
    r_tot = r_n[0] + r_n[1];

    /* 1: When does an event happen? */
    dt = rexp(1 / r_tot);
    t = t + dt;

    /* This is wrong... */
    if ( t > max_t ) {
      dt = dt - (t - max_t);
      t = max_t;
      for ( i = 0; i < n_taxa; i++ )
	len[lineages[i]] += dt;
      break;
    }

    for ( i = 0; i < n_taxa; i++ )
      len[lineages[i]] += dt;

    /* 2: What state does the event happen to? */
    state = unif_rand() > (r_n[0] / r_tot);
    
    /* 3: Given this, what lineage did the event happen to? */
    k = SampleOne(n_i[state]);

    /* Pick a lineage for that state */
    lineage = -1;
    for ( i = 0, j = 0; i < n; i++ ) {
      if ( states[lineages[i]] == state && (j++) == k ) {
	lineage = lineages[i];
	break;
      }
    }
    if ( lineage < 0 )
      error("Something terrible might happen here.");

    /* And pick a type of event */
    /* TODO: This copy is only necessary when directly using
       ProbSampleOne here */
    for ( i = 0; i < 3; i++ ) p[i] = pars[state * 3 + i];
    type = ProbSampleOne(3, p, perm);

    /* Summarise the event */
    if ( type == 0 ) { /* Speciate */
      if ( n_taxa == max_taxa )
	break;

      split[lineage] = 1;
      for ( i = n_entries; i < n_entries + 2; i++ ) {
	/* idx[i] = lineage; */ /* TODO: Don't do this? */
	parent[i] = lineage;
	states[i] = state;
	extinct[i] = 0;
	split[i] = 0;
	start[i] = t;
	len[i] = 0.0;
      }
      (n_i[state])++;
      n_taxa++;
      n_entries += 2;
      /* This can probably be done more efficiently by finding the
	 position of lineage within lineages - which we already know -
	 and pulling everything above that one step down, then
	 appending the two new ones - I may work on this in a
	 bit... */
      for ( i = 0, j = 0; i < n_entries; i++ )
	if ( !extinct[i] && !split[i] )
	  lineages[j++] = i;
    } else if ( type == 1 ) { /* Extinct */
      extinct[lineage] = 1;
      /* See note on lineages in speciation above */
      for ( i = 0, j = 0; i < n_entries; i++ )
	if ( !extinct[i] && !split[i] )
	  lineages[j++] = i;
      (n_i[state]--);
      n_taxa--;
    } else { /* Character change */
      states[lineage] = !state;
      n_i[state]--;
      n_i[!state]++;
      /* This is where I would append history information */
      hist[n_hist][0] = lineage;
      hist[n_hist][1] = state;
      hist[n_hist][2] = !state;
      hist_t[n_hist] = t;
      n_hist++;
    }
  }

  *t_start = t;
  return n_entries;
}

void r_simulate_bisse(double *pars, int *max_taxa, double *max_t, 
		      int *parent, int *states, int *extinct, 
		      int *split, double *start, double *len,
		      int hist[][3], double *hist_t,
		      int *n_entries, int *n, double *t_start) {
  GetRNGstate();
  *n_entries=
    simulate_bisse(pars, *max_taxa, *max_t, parent, states, extinct, 
		   split, start, len, hist, hist_t,
		   *n_entries, *n, t_start);
  PutRNGstate();
}

/*
SEXP simulate_bisse2(SEXP r_pars, SEXP r_max_taxa, SEXP r_max_t, 
		     SEXP r_x0, SEXP r_n) {
  int n = INTEGER(r_n)[0];
  SEXP r_parent, r_states, r_extinct, r_split;
  SEXP r_start, r_len, r_hist, r_hist_t;
  SEXP ret;
  double t_start[1];

  PROTECT(parent = allocVector(INTSXP, n));
  PROTECT(states = allocVector(INTSXP, n));
  PROTECT(extinct = allocVector(LGLSXP, n));
  PROTECT(split = allocVector(LGLSXP, n));

  PROTECT(start = allocVector(INTSXP, n));
  PROTECT(len = allocVector(INTSXP, n));
  PROTECT(hist = allocVector(REALSXP, 3*n));
  PROTECT(hist_t = allocVector(REALSXP, n));

  n_entries = 1;
  t_start[0] = 0.0;

  n = INTEGER(r_n)[0];
  for ( ... ) {
    n_entries = simulate_bisse(REAL(pars),
			       INTEGER(max_taxa)[0],
			       NUMERIC(max_t)[0],
			       INTEGER(parent),
			       INTEGER(states),
			       INTEGER(extinct),
			       INTEGER(split),
			       NUMERIC(start),
			       NUMERIC(len),
			       INTEGER(hist), 
			       NUMERIC(hist_t),
			       n_entries,
			       n,
			       t_start);

    if ( t_start[0] < 0 ) {
      error("need to do this still...");
    }
  }

  PROTECT(ret = allocVector(VECSXP, 10));
  SET_VECTOR_ELT(ret, 1, len);
  SET_VECTOR_ELT(ret, 2, parent);
  SET_VECTOR_ELT(ret, 3, 

  UNPROTECT(8);
  
}
*/
