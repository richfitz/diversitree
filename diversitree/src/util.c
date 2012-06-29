#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>

/* Simplified matrix multiplication, assuming straightforward sizes
   and zeroing the input.  GEMM does:
     Z = alpha X Y + beta Z

   With 
     beta = 0, this gives a fresh calculation of X Y, and with 
     beta = 1, this updates Z to give Z = Z + X Y
*/
void do_gemm(double *x, int nrx, int ncx,
             double *y, int nry, int ncy,
             double *z) {
  const char *trans = "N";
  double alpha = 1.0, beta = 0.0;
  F77_CALL(dgemm)(trans, trans, &nrx, &ncy, &ncx, &alpha,
                  x, &nrx, y, &nry, &beta, z, &nrx);
}

void do_gemm2(double *x, int nrx, int ncx,
	      double *y, int nry, int ncy,
	      double *z) {
  const char *trans = "N";
  double alpha = 1.0, beta = 1.0;
  F77_CALL(dgemm)(trans, trans, &nrx, &ncy, &ncx, &alpha,
                  x, &nrx, y, &nry, &beta, z, &nrx);
}

void do_gemm3(double *x, int nrx, int ncx,
	      double *y, int nry, int ncy,
	      double *z, double beta) {
  const char *trans = "N";
  double alpha = 1.0;
  F77_CALL(dgemm)(trans, trans, &nrx, &ncy, &ncx, &alpha,
                  x, &nrx, y, &nry, &beta, z, &nrx);
}

void r_gemm(double *x, int *nrx, int *ncx,
            double *y, int *nry, int *ncy,
            double *z) {
  do_gemm(x, *nrx, *ncx, y, *nry, *ncy, z);
}

void r_gemm2(double *x, int *nrx, int *ncx,
	     double *y, int *nry, int *ncy,
	     double *z) {
  do_gemm2(x, *nrx, *ncx, y, *nry, *ncy, z);
}

SEXP matrix_to_list(SEXP r_m) {
  SEXP ret, tmp;
  int i, j, k, nr = nrows(r_m), nc = ncols(r_m);
  double *in, *out;

  in = REAL(r_m);

  PROTECT(ret = allocVector(VECSXP, nr));

  for ( i = 0; i < nr; i++ ) {
    /*
      I believe that I don't have to protect agressively here;
      otherwise something like below would be needed.

      PROTECT(tmp = allocVector(REALSXP, nc));
      SET_VECTOR_ELT(ret, i, tmp);
      UNPROTECT(1);
      out = REAL(tmp);

      another option, which definitely does not need garbage
      collection, is:

      SET_VECTOR_ELT(ret, i, allocVector(REALSXP, nc));
      out = REAL(VECTOR_ELT(ret, i));

      which falls somewhere between the two approaches in speed.
      
      However, I've run this under gctorture, and it seems not to
      crash, which is a good sign.
    */
    SET_VECTOR_ELT(ret, i, tmp = allocVector(REALSXP, nc));
    out = REAL(tmp);

    for ( j = 0, k = i; j < nc; j++, k+= nr )
      out[j] = in[k];
  }

  UNPROTECT(1);
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

SEXP getListElementIfThere(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  int i; 
  for ( i = 0; i < length(list); i++ )
    if ( strcmp(CHAR(STRING_ELT(names, i)), str) == 0 ) { 
      elmt = VECTOR_ELT(list, i); 
      break; 
    }

  return elmt;
} 

/* Tree utilities */
int descendants(int node, int *edge, int nedge, int ntip, int *desc);
void descendants_flag(int node, int *edge, int nedge, int ntip, 
		      int *flag);

SEXP r_descendants(SEXP node, SEXP edge, SEXP ntip) {
  int nedge = nrows(edge), *desc = (int *)R_alloc(nedge, sizeof(int));
  int n, *ret_c, node_c = INTEGER(node)[0];
  SEXP ret;
  n = descendants(node_c, INTEGER(edge), nedge,
		  INTEGER(ntip)[0], desc);
  PROTECT(ret = allocVector(INTSXP, n+1));
  ret_c = INTEGER(ret);
  ret_c[0] = node_c;
  memcpy(ret_c + 1, desc, n*sizeof(int));
  UNPROTECT(1);
  return ret;
}

/* if column 1 was sorted, this would be faster... */
int descendants(int node, int *edge, int nedge, int ntip, int *desc) {
  const int *from = edge, *to = edge + nedge;
  int i, n = 0, ni;
  for ( i = 0; i < nedge; i++ ) {
    if ( from[i] == node ) {
      *desc = to[i];
      if ( to[i] > ntip ) /* R indexing... */
	ni = descendants(to[i], edge, nedge, ntip, desc+1) + 1;
      else 
	ni = 1;
      n += ni;
      desc += ni;
    }
  }
  return n;
}

SEXP r_descendants_flag(SEXP node, SEXP edge, SEXP ntip) {
  int nedge = nrows(edge);
  int i, *ret_c, node_c = INTEGER(node)[0];
  int *to = INTEGER(edge) + nedge;
  SEXP ret;
  PROTECT(ret = allocVector(LGLSXP, nedge));
  ret_c = INTEGER(ret);
  for ( i = 0; i < nedge; i++ )
    ret_c[i] = to[i] == node_c;
  descendants_flag(node_c, INTEGER(edge), nedge, INTEGER(ntip)[0],
		   ret_c);
  UNPROTECT(1);
  return ret;
}

void descendants_flag(int node, int *edge, int nedge, int ntip, 
		      int *flag) {
  const int *from = edge, *to = edge + nedge;
  int i;
  for ( i = 0; i < nedge; i++ ) {
    if ( from[i] == node ) {
      flag[i] = 1;
      if ( to[i] > ntip ) /* R indexing... */
	descendants_flag(to[i], edge, nedge, ntip, flag);
    }
  }
}

SEXP r_descendants_idx(SEXP node, SEXP edge, SEXP ntip) {
  int nedge = nrows(edge);
  SEXP ret, flag;
  int *flag_c, *tmp = (int*)R_alloc(nedge, sizeof(int));
  int i, n=0;
  PROTECT(flag = r_descendants_flag(node, edge, ntip));
  flag_c = INTEGER(flag);

  for ( i = 0 ; i < nedge; i++ )
    if ( flag_c[i] )
      tmp[n++] = i + 1;

  PROTECT(ret = allocVector(INTSXP, n));
  memcpy(INTEGER(ret), tmp, n*sizeof(int));
  UNPROTECT(2);

  return ret;
}

SEXP check_ptr_not_null(SEXP extPtr) {
  if ( TYPEOF(extPtr) != EXTPTRSXP )
    error("Recieved non-pointer");
  if ( R_ExternalPtrAddr(extPtr) == NULL )
    error("Recieved NULL pointer");
  return ScalarLogical(1);
}
