
typedef struct {
  int nd;      /* the number of dimensions */
  long nx;     /* number of x positions */
  long ny;     /* floor(nx/2) + 1 */
  double *x;
  complex *y;
  fftw_plan plan_f;
  fftw_plan plan_b;
  int dir; /* DIR_ROWS | DIR_COLS */
} rfftw_plan_real;

/*
  These define the directions that the matrix is organised in.  

  If DIR_COLS is used, then I assume that in the matrix 'x', the
  first 'nx' elements correspond to the first dimension, and the
  second 'nx' elements the second, etc.  This is the way that R stores
  matrices, if you consider the columns.

  If DIR_ROWS is used, then I assume that the matrix 'x' is stored in
  the transpose; the first 'nd' elements correspond to the first
  element of each of the 'nd' dimensions.
 */
#define DIR_COLS 1
#define DIR_ROWS 2

rfftw_plan_real* make_rfftw_plan_real(int nd, int nx, int dir,
				      double *x, complex *y, int flags);
SEXP r_make_rfftw_plan_real(SEXP r_nd, SEXP r_nx, SEXP r_dir);
SEXP r_rfftw_forw(SEXP extPtr, SEXP r_x_in);
SEXP r_rfftw_back(SEXP extPtr, SEXP r_y_in);
SEXP r_get_wisdom();
SEXP r_set_wisdom(SEXP r_wisdom);
