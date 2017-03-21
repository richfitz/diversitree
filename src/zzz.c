#include "zzz.h"
#include <R_ext/Rdynload.h>

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

void R_init_diversitree(DllInfo *dll) {
  R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
  set_sane_gsl_error_handling();
}
