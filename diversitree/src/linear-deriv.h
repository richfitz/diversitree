void build_P_r(int nd, int np, int nt, int with_gr, double *t, 
	       double *d, double *Amat, double *Ainv, 
	       double *Q_d, double *wrk, double *P, double *P_d);
void build_P_c(int nd, int np, int nt, int with_gr, double *t, 
	       Rcomplex *d, Rcomplex *Amat, Rcomplex *Ainv,
	       double *Q_d, Rcomplex *wrk, double *P, double *P_d);
void build_P_1_r(int nd, int np, int with_gr, double t, 
		 double *d, double *Amat, double *Ainv, 
		 double *G, double *wrk,
		 double *P, double *P_d);
void build_P_1_c(int nd, int np, int with_gr, double t,
		 Rcomplex *d, Rcomplex *Amat, Rcomplex *Ainv,
		 Rcomplex *G, Rcomplex *wrk,
		 double *P, double *P_d);
void linear_deriv_tips(int nd, int np, int ntip, int with_gr,
		       int *state, double *P, double *P_d,
		       double *base, double *base_d, double *lq);
void linear_deriv_ints(int nd, int np, int nint, int with_gr,
		       int *order, int *children, 
		       double *P, double *P_d,
		       double *base, double *base_d, double *lq, 
		       double *root, double *root_d);
void linear_deriv_initial_conditions(int nd, int np, int with_gr,
				     double *D1, double *D2,
				     double *D1_d, double *D2_d,
				     double *Dout, double *Dout_d);
double linear_deriv_root(int nd, int np, int nlq, int with_gr,
			 double *root, double *root_d,
			 double *lq, double *gr);
