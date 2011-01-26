/* 
   In this more simple version of the integrator, there is just one
   structure: 
     "quasse_fft"
   This holds everything needed to perform an integration.  It will
   probably be the case that we will want to make many of these to do
   things efficiently, but the C code is agnostic about this, and will
   just provide a few functions to get data in and out of these, as
   well as perform the integrations.

   This version allows there to be multiple nested plans.
*/
typedef struct {
  int n_fft; /* number of different plans */

  int nx;    /* x-extent */
  double dx; /* distance between x's */

  int *nd;    /* number of dimensions for each plan */

  /* Data */
  double *x;
  fftw_complex *y;

  /* Diversification information */
  double *lambda;
  double *mu;

  /* Drift and diffusion parameters */
  double drift;
  double diffusion;

  /* Scratch space */
  double *z; /* Generally stores exp(-rt) */
  double *wrk;

  /* Transform for x propagation */
  rfftw_plan_real **fft;

  /* Kernel information (space propagation) */
  double ny;      /* Fourier space extent */
  fftw_complex *fkern; /* Gaussian kernel transformation */
  int nkl;        /* Kernel width to the left */
  int nkr;        /* Kernel width to the right */
  int npad;       /* nkl + nkr + 1 */
  int ndat;       /* nx - npad */

  /* The kernel itself */
  double  *kern_x;
  fftw_complex *kern_y;
  rfftw_plan_real *kernel;
} quasse_fft;

quasse_fft* make_quasse_fft(int n_fft, int nx, double dx, int *nd, 
			    int flags);
void qf_copy_x(quasse_fft *obj, double *x, int nd, int copy_in);
void qf_copy_ED(quasse_fft *obj, double *x, int idx);
void qf_setup_kern(quasse_fft *obj, double drift, double diffusion,
		   double dt, int nkl, int nkr);
void integrate(quasse_fft *obj, double *lambda, double *mu, 
	       double drift, double diffusion, double *padding);
void propagate_t(quasse_fft *obj, int idx);
void propagate_x(quasse_fft *obj, int idx);
void convolve(rfftw_plan_real *obj, fftw_complex *fy);
int lookup(int x, int *v, int len);
void do_integrate(quasse_fft *obj, int nt, int idx);
