/*
 * TODO: There is still room to do the usual "unique time lengths"
 * optimisation here, but it is not added.  Changes would go in
 * whatever calls build_P() and in linear_deriv_tips(), assuming that
 * only tips are checked for duplicate times.
 */

#include <R.h>
#include <Rinternals.h>
#include "linear-deriv.h"
#include "util-complex.h"
#include "util-matrix.h"

/* This is the main interface function, so it goes up here.

   Arguments:
     r_nd: number of variables (= dimension of Q)
     r_np: number of parameters (for mkn = nd*(nd-1))

     d: vector of eigenvalues; possibly complex.
     Amat: matrix of eigenvectors.  The i'th group within Amat is the
       eigenvector associated with d[i].  Complex iff d is complex.
     Ainv: inverse of Amat.  Complex iff d is complex.

     Q_d: This is np * nd * nd long; the ith set of nd*nd elements is
       dQ/dp_i for the ith parameter.
     t: vector of times (all branch lengths).  NA values are skipped
       over 

     state: vector of tip states.  Numbered so that 0 is the first
       state, and so on.
     order, children: See cache.


*/
SEXP linear_deriv(SEXP r_nd, SEXP r_np, 
		  SEXP d, SEXP Amat, SEXP Ainv,
		  SEXP Q_d, SEXP t,
		  SEXP state, SEXP order, SEXP children,
		  SEXP r_with_gr) {
  SEXP lnlik, gr;
  int nd=INTEGER(r_nd)[0], np=INTEGER(r_np)[0], nt=LENGTH(t);
  int ntip=(nt + 1)/2, nint=ntip - 2;
  int nd2=nd*nd, nd2p=nd2p;
  int with_gr = LOGICAL(r_with_gr)[0];
  int i;

  double *P   = (double*)R_alloc(nt * nd2, sizeof(double));
  double *P_d = (double*)R_alloc(nt * nd2 * np, sizeof(double));

  double *base   = (double*)R_alloc(nt * nd, sizeof(double));
  double *base_d = (double*)R_alloc(nt * nd * np, sizeof(double));

  double *root   = (double*)R_alloc(nd, sizeof(double));
  double *root_d = (double*)R_alloc(nd * np, sizeof(double));
  double *lq     = (double*)R_alloc(nt, sizeof(double));

  /* We may need a real or complex workspace, depending on whether or
     not the eigenvectors are real */
  double *wrk_r;
  Rcomplex *wrk_c;

  /* The lq vector must be zeroed, as it will be referenced but one
     element is never written to... */
  for ( i = 0; i < nt; i++ )
    lq[i] = 0.0;

  if ( isComplex(d) ) {
    wrk_c = (Rcomplex*)R_alloc(nd2p + nd*(1+3*nd) + nd2*(1+np), 
			       sizeof(Rcomplex));
    build_P_c(nd, np, nt, with_gr, REAL(t),
	      COMPLEX(d), COMPLEX(Amat), COMPLEX(Ainv), REAL(Q_d),
	      wrk_c, P, P_d);
  } else {
    wrk_r = (double*)R_alloc(nd2p + nd*(1+3*nd), 
			     sizeof(double));
    build_P_r(nd, np, nt, with_gr, REAL(t),
	      REAL(d), REAL(Amat), REAL(Ainv), REAL(Q_d), 
	      wrk_r, P, P_d);
  }

  linear_deriv_tips(nd, np, ntip, with_gr, INTEGER(state), P, P_d,
		    base, base_d, lq);
  
  linear_deriv_ints(nd, np, nint, with_gr,
		    INTEGER(order), INTEGER(children),
		    P, P_d, base, base_d, lq, root, root_d);

  PROTECT(lnlik = allocVector(REALSXP, 1));

  if ( with_gr ) {
    PROTECT(gr  = allocVector(REALSXP, np));
    REAL(lnlik)[0] = linear_deriv_root(nd, np, nt, with_gr, root, root_d, 
				       lq, REAL(gr));
    setAttrib(lnlik, install("gr"), gr);
  } else {
    REAL(lnlik)[0] = linear_deriv_root(nd, np, nt, with_gr, root, root_d, 
				       lq, NULL);
  }

  UNPROTECT(1 + with_gr);

  return lnlik;
}

SEXP linear_deriv_multi(SEXP r_nd, SEXP r_np, 
			SEXP d, SEXP Amat, SEXP Ainv,
			SEXP Q_d, SEXP t,
			SEXP r_states, SEXP r_w,
			SEXP order, SEXP children,
			SEXP r_with_gr) {
  SEXP r_lnlik, r_gr;
  int nd=INTEGER(r_nd)[0], np=INTEGER(r_np)[0], nt=LENGTH(t);
  int ns=LENGTH(r_w);
  int ntip=(nt + 1)/2, nint=ntip - 2;
  int nd2=nd*nd, nd2p=nd2p;
  int with_gr = LOGICAL(r_with_gr)[0];
  int i, j;
  int *states = INTEGER(r_states), *w = INTEGER(r_w);

  double *P   = (double*)R_alloc(nt * nd2, sizeof(double));
  double *P_d = (double*)R_alloc(nt * nd2 * np, sizeof(double));

  double *base   = (double*)R_alloc(nt * nd, sizeof(double));
  double *base_d = (double*)R_alloc(nt * nd * np, sizeof(double));

  double *root   = (double*)R_alloc(nd, sizeof(double));
  double *root_d = (double*)R_alloc(nd * np, sizeof(double));
  double *lq     = (double*)R_alloc(nt, sizeof(double));

  double lnlik, *gr, *gr_tmp=NULL;

  /* We may need a real or complex workspace, depending on whether or
     not the eigenvectors are real */
  double *wrk_r;
  Rcomplex *wrk_c;

  /* The lq vector must be zeroed, as it will be referenced but one
     element is never written to... */
  for ( i = 0; i < nt; i++ )
    lq[i] = 0.0;

  if ( isComplex(d) ) {
    wrk_c = (Rcomplex*)R_alloc(nd2p + nd*(1+3*nd) + nd2*(1+np), 
			       sizeof(Rcomplex));
    build_P_c(nd, np, nt, with_gr, REAL(t),
	      COMPLEX(d), COMPLEX(Amat), COMPLEX(Ainv), REAL(Q_d),
	      wrk_c, P, P_d);
  } else {
    wrk_r = (double*)R_alloc(nd2p + nd*(1+3*nd), 
			     sizeof(double));
    build_P_r(nd, np, nt, with_gr, REAL(t),
	      REAL(d), REAL(Amat), REAL(Ainv), REAL(Q_d), 
	      wrk_r, P, P_d);
  }

  PROTECT(r_lnlik = allocVector(REALSXP, 1));
  PROTECT(r_gr  = allocVector(REALSXP, np));
  lnlik = 0;
  gr = REAL(r_gr);
  if ( with_gr ) {
    gr_tmp = (double*)R_alloc(np, sizeof(double));
    for ( j = 0; j < np; j++ )
      gr[j] = 0.0;
  }

  /* Here, we reuse the P and P_d matrices to compute the individual
     likelihoods */
  for ( i = 0; i < ns; i++ ) {
    linear_deriv_tips(nd, np, ntip, with_gr, states + i*ntip, P, P_d,
		      base, base_d, lq);

    linear_deriv_ints(nd, np, nint, with_gr,
		      INTEGER(order), INTEGER(children),
		      P, P_d, base, base_d, lq, root, root_d);

    lnlik += w[i] * linear_deriv_root(nd, np, nt, with_gr, root, root_d, 
				      lq, gr_tmp);

    if ( with_gr )
      for ( j = 0; j < np; j++ )
	gr[j] += w[i] * gr_tmp[j];
  }

  REAL(r_lnlik)[0] = lnlik;
  if ( with_gr )
    setAttrib(r_lnlik, install("gr"), r_gr);

  UNPROTECT(2);

  return r_lnlik;
}

/* Given the eigensystem and vector of times, compute the transition
   probability density matrix, P = Pij, for each branch.  With the
   help of Q_d = d(Q)/dp, compute d(P)/dp for each branch.

   The workspace vector 'wrk' here needs to be of length
     nd*nd*np + nd(1+3nd)
   assigned for
     G: nd*nd*np
     wrk_1: nd(1+3*nd)
 */
void build_P_r(int nd, int np, int nt, int with_gr, double *t, 
	       double *d, double *Amat, double *Ainv, 
	       double *Q_d, double *wrk, double *P, double *P_d) {
  const int nd2 = nd * nd, nd2p = nd * nd * np;
  int i, off_P=0, off_P_d=0;
  double *G = wrk, *wrk_1 = wrk + nd2p;

  /* Compute the G matrix out of Q_d */
  if( with_gr )
    for ( i = 0; i < np; i++ )
      mult_mmm(nd, Ainv, Q_d + i * nd2, Amat, 
	       G + i * nd2, wrk_1);

  for ( i = 0; i < nt; i++, off_P += nd2, off_P_d += nd2p )
    build_P_1_r(nd, np, with_gr, t[i], d, Amat, Ainv, G, wrk_1, 
		P + off_P, P_d + off_P_d);
}

/* This is just like build_P_r(), but we have complex values for d,
   Amat, and Ainv.  The complex part will be stripped off all results
   so that P and P_d are real.

   The workspace vector here is longer than the real version:
     nd*nd*np + nd(1+3nd) + nd*nd(1+np)
   assigned for
     G: nd*nd*np
     wrk_1: nd(1+3nd) + nd*nd(1+np) (see documentation for
   build_P_1_c())
*/
void build_P_c(int nd, int np, int nt, int with_gr, double *t, 
	       Rcomplex *d, Rcomplex *Amat, Rcomplex *Ainv,
	       double *Q_d, Rcomplex *wrk, double *P, double *P_d) {
  const int nd2 = nd * nd, nd2p = nd * nd * np;
  int i, off_P=0, off_P_d=0;
  Rcomplex *G = wrk, *wrk_1 = wrk + nd2p;
  Rcomplex *zQ_d = wrk_1 + nd2;

  /* First up, we have to copy Q_d into a complex vector */
  for ( i = 0; i < nd2p; i++ ) {
    zQ_d[i].r = Q_d[i];
    zQ_d[i].i = 0.0;
  }

  /* Compute the G matrix out of Q_d */
  if( with_gr )
    for ( i = 0; i < np; i++ )
      zmult_mmm(nd, Ainv, zQ_d + i * nd2, Amat, 
		G + i * nd2, wrk_1);

  /* Then run all the usual eigenvalue stuff, but through a complex
     version. */
  for ( i = 0; i < nt; i++, off_P += nd2, off_P_d += nd2p )
    build_P_1_c(nd, np, with_gr, t[i], d, Amat, Ainv, G, wrk_1, 
		P + off_P, P_d + off_P_d);
}

/*
 * Inputs:
 *   nd: number of variables (dimension of Q)
 *   np: number of parameters (for mkn, nd*(nd - 1)
 *   t: time length for the calculation
 *   d: vector of eigenvalues
 *   Amat: matrix of eigenvectors
 *   Ainv: inverse of Amat
 *   G: (tricky to get right) matrix of derivatives of each Q matrix
 *      with respect to each parameter.  This is np * nd * nd long;
 *      each set of nd*nd elements is dQ/dp_i for a parameter p_i.
 *   wrk: some scratch space.  Currently, I need
 *      - nd for exp(d*t)
 *      - nd*nd for V
 *      - nd*nd for V * G
 *      - nd*nd for mult_mmm() scratch
 *      for a total of nd(1 + 3*nd)
 * Outputs:
 *   P: Pij = Amat exp(diag(d)t) Amat^{-1}
 *   P_d: same format as G; each nd^2 elements is dpij/dp_k for a
 *     parameter p_k
 */
void build_P_1_r(int nd, int np, int with_gr, double t, 
		 double *d, double *Amat, double *Ainv, 
		 double *G, double *wrk,
		 double *P, double *P_d) {
  const int nd2 = nd * nd;
  int i, j, k1, k2;
  double *edt = wrk + nd2;
  double *V   = edt + nd;
  double *VG  = V + nd2;
  double *Gi = G, *P_di = P_d;

  if ( ISNA(t) ) return;

  for ( i = 0; i < nd; i++ )
    edt[i] = exp(d[i] * t);

  mult_mdm(nd, Amat, edt, Ainv, P, wrk);

  if ( with_gr ) {
    for ( i = 0; i < nd; i++ ) {
      k1 = k2 = i * nd + i;
      V[k1] = t * edt[i];
      for ( j = i + 1; j < nd; j++ ) {
	k1++;
	k2+=nd;
	V[k1] = V[k2] = (edt[i] - edt[j]) / (d[i] - d[j]);
      }
    }

    for ( i = 0; i < np; i++, Gi += nd2, P_di += nd2 ) {
      for ( j = 0; j < nd2; j++ )
	VG[j] = V[j] * Gi[j];
      mult_mmm(nd, Amat, VG, Ainv, P_di, wrk);
    }
  }
}

/* vs the normal build_P_1_r, here are the changes:
   d, Amat, Ainv are complex (duh)
   G is complex (just made it so above)
   
   at least *some* of wrk is complex

   wrk:
     nd for exp(d*t) [complex]
     nd * nd for V [complex?]
     nd * nd for VG [complex]
     nd * nd for zmult_mmm() scratch (complex)
   so we have nd(1+3nd) as above. PLUS:
     nd * nd for tmp_P
     nd * nd * np for tmp_P_d
   for a total of nd(1+3nd) + nd*nd(1+np)
 */
void build_P_1_c(int nd, int np, int with_gr, double t, 
		 Rcomplex *d, Rcomplex *Amat, Rcomplex *Ainv,
		 Rcomplex *G, Rcomplex *wrk,
		 double *P, double *P_d) {
  const int nd2 = nd * nd;
  int i, j, k1, k2;
  Rcomplex *edt = wrk + nd2;
  Rcomplex *V   = edt + nd;
  Rcomplex *VG  = V + nd2;
  Rcomplex *tmp_P = VG + nd2;
  Rcomplex *tmp_P_d = tmp_P + nd2;
  Rcomplex *Gi = G;
  double *P_di = P_d;

  if ( ISNA(t) ) return;

  for ( i = 0; i < nd; i++ )
    edt[i] = z_exp(z_scal(d[i], t));

  zmult_mdm(nd, Amat, edt, Ainv, tmp_P, wrk);
  for ( i = 0; i < nd2; i++ )
    P[i] = tmp_P[i].r;

  if ( with_gr ) {
    for ( i = 0; i < nd; i++ ) {
      k1 = k2 = i * nd + i;
      V[k1] = z_scal(edt[i], t);
      for ( j = i + 1; j < nd; j++ ) {
	k1++;
	k2+=nd;
	V[k1] = V[k2] = z_divide(z_minus(edt[i], edt[j]),
				 z_minus(d[i], d[j]));
      }
    }

    for ( i = 0; i < np; i++, Gi += nd2, P_di += nd2 ) {
      for ( j = 0; j < nd2; j++ )
	VG[j] = z_times(V[j], Gi[j]);

      zmult_mmm(nd, Amat, VG, Ainv, tmp_P_d, wrk);

      for ( j = 0; j < nd2; j++ )
	P_di[j] = tmp_P_d[j].r;
    }
  }
}

/* Given P and P_d arrays, compute base probabilities and gradients.
   This is done separately from the internal nodes, as when all states
   are known this is really easy 

   TODO: For unknown states (indicated perhaps with state[i] < 0?), we
   must multiply Pij * {1}/nd, which is the same as summing over the
   rows.  This might be better to do with mult_mv(), but it should
   still be easy enough to do.
*/
void linear_deriv_tips(int nd, int np, int ntip, int with_gr,
		       int *state, double *P, double *P_d,
		       double *base, double *base_d, double *lq) {
  const int nd2 = nd * nd, nd2p = nd * nd * np;
  int i, j, k;
  double *P_i, *P_d_i, *base_i, *base_d_i;
  double q;

  for ( i = 0; i < ntip; i++ ) {
    q = 0;
    P_i = P + i * nd2 + state[i] * nd;
    base_i = base + i * nd;

    for ( j = 0; j < nd; j++ )
      q += P_i[j];

    for ( j = 0; j < nd; j++ )
      base_i[j] = P_i[j] / q;

    lq[i] = log(q);

    /* At this point, P_d, offset to a particular branch base
       contains np blocks of nd*nd elements.  The kth of these blocks
       is the derivative of Pij with respect to the kth parameter.  As
       above, we want the state[i]'th column from each of these blocks
       to form the columns in base_d_i 

       Leaving here, base_d has structure so that each 'nd*np'
       elements are the derivatives for one branch.  Within this, each
       nd elements are the derivatives of all variables with respect
       to a particular parameter.
    */
    if ( with_gr ) {
      P_d_i = P_d + i * nd2p + state[i] * nd;
      base_d_i = base_d + i * nd * np;

      for ( k = 0; k < np; k++ ) {
	for ( j = 0; j < nd; j++ ) 
	  base_d_i[j] = P_d_i[j] / q;
	P_d_i += nd2;
	base_d_i += nd;
      }
    }
  }
}

/* Internal nodes */
void linear_deriv_ints(int nd, int np, int nint, int with_gr,
		       int *order, int *children, 
		       double *P, double *P_d,
		       double *base, double *base_d, double *lq, 
		       double *root, double *root_d) {
  const int ndp = nd * np, nd2 = nd * nd, nd2p = nd * nd * np;
  int i, j, idx, idx1, idx2;
  double *init   = root;
  double *init_d = root_d;
  double *P_i, *P_d_i, *base_i, *base_d_i, q;

  for ( i = 0; i < nint; i++ ) {
    idx = order[i];               /* self         */
    idx1 = children[idx * 2];     /* first child  */
    idx2 = children[idx * 2 + 1]; /* second child */
    
    P_i    = P    + idx * nd2;
    base_i = base + idx * nd;

    linear_deriv_initial_conditions(nd, np, with_gr,
				    base + idx1 * nd,
				    base + idx2 * nd,
				    base_d + idx1 * ndp,
				    base_d + idx2 * ndp,
				    init, init_d);
    
    mult_mv(nd, P_i, init, 0.0, base_i);

    q = 0;
    for ( j = 0; j < nd; j++ )
      q += base_i[j];
    for ( j = 0; j < nd; j++ )
      base_i[j] /= q;
    lq[idx] = log(q);

    /* Here, we want to compute
         d(base)   dP               d(init)
         ------ = ---- * init + P * ------
          dp_j    dp_j               dp_j
       where "*" is matrix multiplication, and the second element in
       each of these is a vector */
    if ( with_gr ) {
      P_d_i    = P_d    + idx * nd2p;
      base_d_i = base_d + idx * ndp;

      for ( j = 0; j < np; j++ ) {
	mult_mv(nd, P_d_i + j * nd2, init, 0.0, base_d_i + j * nd);
	mult_mv(nd, P_i, init_d + j * nd,  1.0, base_d_i + j * nd);
      }
      for ( j = 0; j < ndp; j++ )
	base_d_i[j] /= q;
    }
  }

  /* Root node */
  idx = order[nint];            /* self         */
  idx1 = children[idx * 2];     /* first child  */
  idx2 = children[idx * 2 + 1]; /* second child */
  linear_deriv_initial_conditions(nd, np, with_gr,
				  base + idx1 * nd,
				  base + idx2 * nd,
				  base_d + idx1 * ndp,
				  base_d + idx2 * ndp,
				  root, root_d);
}

/* Compute the initial conditions at a node

   This is the place where other models would probably deviate; here I
   am assuming that the node calculation is DL * DR for the vectors DL
   and DR (which gives the DL_d * DR + DL * DR_d for the gradient).

   If different calculations occur at the node, then this would need
   to change.
*/
void linear_deriv_initial_conditions(int nd, int np, int with_gr,
				     double *D1, double *D2,
				     double *D1_d, double *D2_d,
				     double *Dout, double *Dout_d) {
  int i, j = 0, k;

  for ( i = 0; i < nd; i++ )
    Dout[i] = D1[i] * D2[i];

  if ( with_gr ) {
    for ( k = 0; k < np; k++ )
      for ( i = 0; i < nd; i++, j++ )
	Dout_d[j] = D1[i] * D2_d[j] + D2[i] * D1_d[j];
  }
}

/* This assumes a flat root, but it would be trivial to extend it to
   any vector of probabilities that are independent of the parameter
   values */
double linear_deriv_root(int nd, int np, int nlq, int with_gr,
			 double *root, double *root_d,
			 double *lq, double *gr) {
  double sum_root=0.0, lnlik, *root_d_k, tmp;
  int i, k;
  for ( i = 0; i < nd; i++ )
    sum_root += root[i];

  lnlik = log(sum_root / nd); /* here is the prior */
  for ( i = 0; i < nlq; i++ )
    lnlik += lq[i];

  if ( with_gr ) {
    root_d_k = root_d;
    for ( k = 0; k < np; k++ ) {
      tmp = 0.0;
      for ( i = 0; i < nd; i++ )
	tmp += root_d_k[i];
      gr[k] = tmp / sum_root;
      root_d_k += nd;
    }
  }
  
  return lnlik;
}
