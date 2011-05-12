RCvodesObj* make_cvodes_fwd(int neq, int np, CVRhsFn rhs, 
			    CVSensRhs1Fn sens1, double rtol,
			    double *atol);
void cvodes_fwd_cleanup(RCvodesObj *obj);
int cvodes_fwd_run(RCvodesObj *obj, double *y0, double *ys0, 
		   double *times, int nt,
		   double *retY, double *retS);
void copy_yS_out(int neq, int nS, int nt, N_Vector *yS, double *out);
void copy_yS_in(int neq, int nS, int nt, N_Vector *yS, double *in);
