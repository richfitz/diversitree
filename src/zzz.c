#include "zzz.h"
#include <R_ext/Rdynload.h>
#include <Rversion.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef R_CallDef[] = {
   CALLDEF(r_all_branches_cont, 2),
   CALLDEF(r_asr_marginal_mkn, 7),
   CALLDEF(r_branches_mkn_expokit, 8),
   CALLDEF(r_check_ptr_not_null, 1),
   CALLDEF(r_descendants, 3),
   CALLDEF(r_descendants_flag, 3),
   CALLDEF(r_descendants_idx, 3),
   CALLDEF(r_do_asr_joint, 7),
   CALLDEF(r_do_integrate, 9),
   CALLDEF(r_do_tips, 9),
   CALLDEF(r_dt_cont_reset_tips, 2),
   CALLDEF(r_get_vals_cont, 1),
   CALLDEF(r_get_wisdom, 0),
   CALLDEF(r_hdr, 3),
   CALLDEF(r_make_dt_obj_cont, 3),
   CALLDEF(r_make_quasse_fft, 4),
   CALLDEF(r_matrix_to_list, 1),
   CALLDEF(r_set_wisdom, 1),
   CALLDEF(r_smkn_alloc, 2),
   CALLDEF(r_smkn_scm_run_all, 7),
   {NULL, NULL, 0}
};

static R_CMethodDef R_CDef[] = {
  {"r_mkn_core",       (DL_FUNC) &r_mkn_core,       8,  NULL, NULL},
  {"r_simulate_bisse", (DL_FUNC) &r_simulate_bisse, 15, NULL, NULL},
  {NULL, NULL, 0, NULL, NULL}
};

// ddexpmv
// dsexpmvi
// dexpmf
static R_FortranMethodDef R_FortranDef[] = {
  {"f_bucexp",  (DL_FUNC) &F77_SUB(bucexp),  14, NULL, NULL},
  {"f_bucexpl", (DL_FUNC) &F77_SUB(bucexpl), 19, NULL, NULL},
  {"f_nucexp",  (DL_FUNC) &F77_SUB(nucexp),  18, NULL, NULL},
  {"f_nucexpl", (DL_FUNC) &F77_SUB(nucexpl), 23, NULL, NULL},
  //
  {"f_ddexpmv",  (DL_FUNC) &F77_SUB(ddexpmv),   6, NULL, NULL},
  {"f_dsexpmvi", (DL_FUNC) &F77_SUB(dsexpmvi), 12, NULL, NULL},
  {"f_dexpmf",   (DL_FUNC) &F77_SUB(dexpmf),    5, NULL, NULL},
  {NULL, NULL, 0, NULL, NULL}
};

void R_init_diversitree(DllInfo *dll) {
  R_registerRoutines(dll, R_CDef, R_CallDef, R_FortranDef, NULL);
#if defined(R_VERSION) && R_VERSION >= R_Version(3, 3, 0)
  //  R_useDynamicSymbols(dll, FALSE);
#endif
  set_sane_gsl_error_handling();
}
