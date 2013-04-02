void mult_mv(int n, double *A, double *v, double beta, double *out);
void mult_mm(int n, double *A, double *B, double *out);
void mult_md(int n, double *A, double *d, double *out);
void mult_mmm(int n, double *A, double *B, double *C, 
	      double *out, double *wrk);
void mult_mdm(int n, double *A, double *d, double *C, double *out,
	      double *wrk);
