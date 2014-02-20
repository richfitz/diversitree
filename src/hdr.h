double RSRC_Brent_fmin(double ax, double bx, double (*f)(double, void *),
		       void *info, double tol);
typedef struct hdr_struct
{
  dt_spline *spline;
  double     w; /* width: 1-alpha */
} hdr_struct, *HDRStruct;
