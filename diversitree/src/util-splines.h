typedef struct {
  int nx;
  double *x;
  double *y;
  double *b;
  double *c;
  double *d;
} dt_spline;

dt_spline* make_dt_spline(int nx, double *x, double *y, int deriv);
void cleanup_dt_spline(dt_spline *obj);
double dt_spline_eval1(dt_spline *obj, double u);
void dt_spline_eval(dt_spline *obj, double *u, int nu, double *ret);
