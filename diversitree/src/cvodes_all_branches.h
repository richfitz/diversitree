/* Typedef for the initial conditions function.

   For what I am doing here, see:

   cvodes.h:285 (function typedef)
   cvodes_impl.h:84 (containing struct def)
   cvodes.c:849 (assignment)
   cvodes.c:2691 (use)
 */
typedef int (*DtIcFun)(int neq, double *vars_l, double *vars_r,
		       double *pars, double t, double *vars_out);

/* I need a good data structure to hold everything */
typedef struct {
  int neq;
  int np;

  RCvodesObj *integrator;
  DtIcFun    ic;

  /* Control underflow compensation */
  int *comp_idx;
  int  comp_n;
  double eps;

  /* Space requirements */
  int n_out;

  /* Tips */
  int      tip_types;  /* number of different tip classes   */
  int     *tip_n;      /* number of species in each class   */
  double **tip_y;      /* initial conditions for each class */
  double **tip_len;    /* branch lengths                    */
  int    **tip_target; /*  target indices (base 0)          */

  /* Internal branches */
  int     n_int;     /* number of internal branches         */
  int    *order;     /* order to process internal branches  */
  int    *children;  /* children (documentation needed)     */
  double *len;       /* vector of branch lengths            */
  double *depth;     /* vector of branch depths             */
  int     root;      /* root node                           */

  /* Storage of results */
  double *init;
  double *base;
  double *lq;
} dt_obj;

SEXP getListElement(SEXP list, const char *str);
void dt_setup_tips(dt_obj *obj, SEXP cache);
void dt_setup_internal(dt_obj *obj, SEXP cache);

/* New version */
void dt_cvodes_run(RCvodesObj *obj, double *y_in,
		   double *len, int nt, double t0,
		   int *comp_idx, int comp_n, double eps,
		   int *target, double *base, double *lq);
double dt_cvodes_run_1(RCvodesObj *obj, double t0, double t1, 
		       int *comp_idx, int comp_n, double eps);
double dt_cvodes_run_multi(RCvodesObj *obj, double *y0, 
			   double t0, double t1, 
			   int *comp_idx, int comp_n, double eps,
			   int depth);

