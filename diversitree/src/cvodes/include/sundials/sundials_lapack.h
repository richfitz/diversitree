#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#define dcopy_f77  F77_CALL(dcopy)
#define dscal_f77  F77_CALL(dscal)
#define dgemv_f77  F77_CALL(dgemv)
#define dtrsv_f77  F77_CALL(dtrsv)
#define dsyrk_f77  F77_CALL(dsyrk)

#define dgbtrf_f77 F77_CALL(dgbtrf)
#define dgbtrs_f77 F77_CALL(dgbtrs)
#define dgetrf_f77 F77_CALL(dgetrf)
#define dgetrs_f77 F77_CALL(dgetrs)
#define dgeqp3_f77 F77_CALL(dgeqp3)
#define dgeqrf_f77 F77_CALL(dgeqrf)
#define dormqr_f77 F77_CALL(dormqr)
#define dpotrf_f77 F77_CALL(dpotrf)
#define dpotrs_f77 F77_CALL(dpotrs)
