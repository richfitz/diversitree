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

  /* Target (purely to deal with Q) */
  int *target;

  int *nonnegative;

  /* Time domain */
  double t_range[2];

  /* Support for q */
  int  idx_q;    /* index of the q matrix, negative if not present  */
  int  k;        /* dimension of the q matrix                       */
  int *const_q;  /* array[k]; const_q[i] == 1 if row i all constant */
  double *q_out; /* pointer to q within p_out                       */
  
  dt_spline *spline_data;

  //void *spline_data; /* spline data */
} dt_time_machine;

/* These are the useful interface functions for other C bits */
void init_time_machine(dt_time_machine *obj, double *pars);
void run_time_machine(dt_time_machine *obj, double t);

/* These are forward references */
void normalise_q(dt_time_machine *obj, int is_const);

/* 1. Here are possible types of functions, and their constants */
#define T_CONSTANT      0
#define T_LINEAR        1
#define T_STEPF         2
#define T_SIGMOID       3
#define T_SPLINE        4
#define T_SPLINE_LINEAR 5
