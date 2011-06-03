/* This is a stripped down version of the MuSSE SCM object. */
typedef struct {
  int k;
  double *pars;

  /* This is a vector containing the rate at which things happen to
     state i */
  double *r;
  double r_tot;

  /* This is a matrix with the cumulative probabilities for each
     state, and the permutation matrix (see SampleOneProb).  Each of
     the k sets of k-1 elements applies to one of the states */
  double *cp;
  int *perm;

  /* This is different to the musse version */
  int state;

  /* Character history */
  int n_hist;
  int n_hist_max;
  int *hist_from;
  int *hist_to;
  double *hist_t;
} smkn_info;

smkn_info* smkn_alloc(int k, int n_hist);
void smkn_cleanup(smkn_info *obj);
void smkn_set_pars(smkn_info* obj, double *pars);

void smkn_set_pars(smkn_info* obj, double *pars);

int smkn_scm_run(smkn_info *obj, double len, 
		 int state_beg, int state_end);
double smkn_sim(smkn_info *obj, int x0, double len);
void smkn_init(smkn_info *obj, int x0);
int smkn_pick_state(smkn_info *obj, int state);
void smkn_evolve(smkn_info* obj, int state, double t, int state_to);
void smkn_grow_hist(smkn_info *obj);

SEXP smkn_slim(SEXP obj);
