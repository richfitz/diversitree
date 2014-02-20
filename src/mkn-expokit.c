#include <R.h>
#include <Rinternals.h>

#include "mkn-expokit.h"

SEXP r_branches_mkn_expokit(SEXP Q, SEXP ia, SEXP ja, SEXP qnorm, 
			    SEXP t, SEXP v, SEXP rm, SEXP tol) {
  int m = INTEGER(rm)[0], n = LENGTH(v);
  const int nmax = 1024, mmax = 30;
  int lwsp = nmax*(mmax+1)+nmax+5*(mmax+2)*(mmax+2)+7,liwsp=nmax+2;
  double * wsp = (double*) R_alloc(lwsp,  sizeof(double));
  int    *iwsp = (int   *) R_alloc(liwsp, sizeof(int));
  int nt = LENGTH(t), nz = LENGTH(ia), itrace = 0, iflag = 0;
  double *cvals, *cscal, *vals_i, tot;
  SEXP ret, vals, scal;
  int i, j;

  if ( m >= n ) {
    warning("Decreasing 'm' to %d", n-1);
    m = n-1;
  }

  PROTECT(ret  = allocVector(VECSXP, 2));
  PROTECT(vals = allocMatrix(REALSXP, n, nt));
  PROTECT(scal = allocVector(REALSXP, nt));
  SET_VECTOR_ELT(ret, 0, scal);
  SET_VECTOR_ELT(ret, 1, vals);
  cvals = REAL(vals);
  cscal = REAL(scal);

  F77_CALL(dgexpvi)(&n, &m,
		    REAL(t), &nt, REAL(v),
		    cvals,
		    REAL(tol), REAL(qnorm),
		    INTEGER(ia), INTEGER(ja), REAL(Q), &nz,
		    wsp, &lwsp, iwsp, &liwsp, &itrace, &iflag);

  if ( iflag != 0 ) {
    if ( iflag == -42 )
      error("expokit failed, but I have no idea why -- try ode instead?");
    else
      error("expokit failed with flag %d\n", iflag);
  }

  for ( i = 0; i < nt; i++ ) {
    tot = 0.0;
    vals_i = cvals + i*n;
    for ( j = 0; j < n; j++ )
      tot += vals_i[j];
    for ( j = 0; j < n; j++ )
      vals_i[j] /= tot;
    cscal[i] = log(tot);
  }

  UNPROTECT(3);
  return ret;
}

