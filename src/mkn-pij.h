void asr_normalise(int n_states, double *vals);
void asr_marginal_mkn_1(int k, int node, int root,
			int *parent, int *children,
			double *pij, double *branch_init, 
			double *branch_base, double *lq);
