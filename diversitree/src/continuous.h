/* The DtIcFun is duplicated from elsewhere, and DtBrFun is new */
typedef int (*DtIcFun)(int neq, double *vars_l, double *vars_r,
		       double *pars, double t, double *vars_out);
typedef double (*DtBrFun)(double *vars_in, double len, double *pars, 
			  double t0, int idx, double *vars_out);

typedef struct {
  int neq;
  int np;

  DtBrFun br;
  DtIcFun ic;

  /* Space requirements */
  int n_out;

  /* Tips */
  /* In theory, two different sorts of tips should be possible here,
     but I'm just doing the simple one */
  int     n_tip;      /* number of tips                     */
  int    *tip_target; /* target within init & base (base 0) */
  double *tip_len;    /* branch lengths of initial tips     */
  double *tip_y;      /* initial conditions for each tip    */

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
} dt_obj_cont;

void dt_cont_setup_tips(dt_obj_cont *obj, SEXP cache);
void dt_cont_setup_internal(dt_obj_cont *obj, SEXP cache);
SEXP r_dt_cont_reset_tips(SEXP extPtr, SEXP tips);
