void do_gemm(double *x, int nrx, int ncx,
             double *y, int nry, int ncy,
             double *z);
void do_gemm2(double *x, int nrx, int ncx,
	      double *y, int nry, int ncy,
	      double *z);
void do_gemm3(double *x, int nrx, int ncx,
	      double *y, int nry, int ncy,
	      double *z, double beta);
SEXP getListElement(SEXP list, const char *str);
SEXP getListElementIfThere(SEXP list, const char *str);
