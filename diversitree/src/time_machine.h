typedef struct {
  /* Number, and array, of high level (hyper) parameters */
  int np_in;
  double *p_in;
  /* Number, and array, of underlying model parameters */
  int np_out;
  double *p_out;
  /* Number of functions, usually the same as np_out */
  int nf;

  /* Array[np_out] of types of functions.  In non-Q models, this will
     be length np_out.  However, for models with Q, this is going to
     change. */
  int *types; /* type of each *output* parameter */

  /* start[i] contains the starting position of function i within
     p_in; we would call f[i] with parameters &p_in[start[i]] */
  int *start;

  //void *spline_data; /* spline data */
} dt_time_machine;

/* These are the useful interface functions for other C bits */
void init_time_machine(dt_time_machine *obj, double *pars);
void run_time_machine(dt_time_machine *obj, double t);
SEXP r_set_tm_bd_t2(SEXP extPtr);
