#include <R.h>
#include <Rinternals.h>

// util.c
SEXP r_check_ptr_not_null(SEXP extPtr);
SEXP r_matrix_to_list(SEXP r_m);
SEXP r_descendants(SEXP node, SEXP edge, SEXP ntip);
SEXP r_descendants_flag(SEXP node, SEXP edge, SEXP ntip);
SEXP r_descendants_idx(SEXP node, SEXP edge, SEXP ntip);
SEXP set_sane_gsl_error_handling();

// hdr.c
SEXP r_hdr(SEXP x, SEXP y, SEXP alpha);

// continuous.c
SEXP r_all_branches_cont(SEXP extPtr, SEXP r_pars);
SEXP r_dt_cont_reset_tips(SEXP extPtr, SEXP tips); // in .h
SEXP r_get_vals_cont(SEXP extPtr);
SEXP r_make_dt_obj_cont(SEXP cache, SEXP r_ic, SEXP r_br);

// mkn-pij
SEXP r_asr_marginal_mkn(SEXP r_k, SEXP r_pars, SEXP r_nodes,
                        SEXP cache, SEXP res,
                        SEXP root_f, SEXP rho);

// mkn_expokit
SEXP r_branches_mkn_expokit(SEXP Q, SEXP ia, SEXP ja, SEXP qnorm,
                            SEXP t, SEXP v, SEXP rm, SEXP tol);

// asr-joint.c
SEXP r_do_asr_joint(SEXP r_k, SEXP r_order, SEXP r_parent,
                    SEXP r_li, SEXP r_pij, SEXP r_root_p, SEXP r_as_01);

// quasse-eqs-fftC
SEXP r_do_integrate(SEXP extPtr, SEXP vars, SEXP lambda, SEXP mu,
                    SEXP drift, SEXP diffusion, SEXP nt, SEXP dt,
                    SEXP padding);
SEXP r_do_tips(SEXP extPtr, SEXP vars, SEXP lambda, SEXP mu,
               SEXP drift, SEXP diffusion, SEXP nt, SEXP dt,
               SEXP padding);
SEXP r_make_quasse_fft(SEXP r_nx, SEXP r_dx, SEXP r_nd, SEXP r_flags);

// rfftw.h
SEXP r_get_wisdom();
SEXP r_set_wisdom(SEXP r_wisdom);

// scm-mkn.c
SEXP r_smkn_alloc(SEXP k, SEXP n_hist);
SEXP r_smkn_scm_run_all(SEXP extPtr, SEXP pars, SEXP r_len,
                        SEXP r_state_beg, SEXP r_state_end,
                        SEXP r_as01, SEXP r_slim);

// mkn-pij.c
void r_mkn_core(int *k, int *n, int *order, int *children, double *pij,
                double *branch_init, double *branch_base, double *lq);

// simulate-bisse.c
void r_simulate_bisse(double *pars, int *max_taxa, double *max_t,
                      int *parent, int *states, int *extinct,
                      int *split, double *start, double *len,
                      int hist[][3], double *hist_t, int *n_info,
                      int *n, double *t_start, int *verbose);

void F77_SUB(ddexpmv)(double *Q, int *n, double *v, double *t,
                      double *out, int *iflag);
void F77_SUB(dsexpmvi)(double *Q, int *n, int *ia, int *ja, int *nz,
                       double *qnorm, double *v, double *t, int *lt,
                       double *tol, double *out, int *iflag);
void F77_SUB(dexpmf)(double *Q, int *n, double *t, double *out, int *iflag);

void F77_SUB(bucexp)(int *nt,
                     double *la0,
                     double *la1,
                     double *mu0,
                     double *mu1,
                     double *q01,
                     double *q10,
                     double *t,
                     int *lt,
                     double *scal,
                     double *tol,
                     int *m,
                     double *w,
                     int *iflag);
void F77_SUB(bucexpl)(int *nt,
                      double *la0,
                      double *la1,
                      double *mu0,
                      double *mu1,
                      double *q01,
                      double *q10,
                      double *t,
                      int *lt,
                      int *ti,
                      int *Nc,
                      int *nsc,
                      int *k,
                      int *lc,
                      double *scal,
                      double *tol,
                      int *m,
                      double *ans,
                      int *iflag);
void F77_SUB(nucexp)(int *nt,
                     double *la0,
                     double *la1,
                     double *mu0,
                     double *mu1,
                     double *q01,
                     double *q10,
                     double *p0c,
                     double *p0a,
                     double *p1c,
                     double *p1a,
                     double *t,
                     int *lt,
                     double *scal,
                     double *tol,
                     int *m,
                     double *w,
                     int *iflag);
void F77_SUB(nucexpl)(int *nt,
                      double *la0,
                      double *la1,
                      double *mu0,
                      double *mu1,
                      double *q01,
                      double *q10,
                      double *p0c,
                      double *p0a,
                      double *p1c,
                      double *p1a,
                      double *t,
                      int *lt,
                      int *ti,
                      int *Nc,
                      int *nsc,
                      int *k,
                      int *lc,
                      double *scal,
                      double *tol,
                      int *m,
                      double *ans,
                      int *iflag);
