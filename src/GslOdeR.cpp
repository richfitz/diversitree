#include "GslOdeR.h"

GslOdeR::GslOdeR(SEXP fun_, SEXP env_, size_t size) : 
  GslOdeBase(size), fun(fun_), env(env_), pars(R_NilValue) {
  // TODO: check that fun is a function and that env is an
  // environment.  Could probably use Rcpp types here to help?
}

void GslOdeR::set_pars(SEXP pars_) {
  pars = pars_;
}

void GslOdeR::clear_pars() {
  pars = R_NilValue;
}

// TODO: Not sure that this actually belongs here, rather than in
// derivs() below.
SEXP GslOdeR::target(double t, SEXP y) {
  return Rf_eval(Rf_lang4(fun, Rf_ScalarReal(t), y, pars), env);
}

void GslOdeR::derivs(double t, const double y[], double dydt[]) {
  // It is possible that we could allocate the space at construction
  // (though I think that allocVector could lead to eventual garbage
  // collection as we can't protect it).  Alternatively,
  // Rcpp::NumericVector followed by Rcpp::wrap might be less ugly.
  // That would require a change in target() to return
  // Rcpp::NumericVector though, and it won't do that neatly.

  // 1. Create R vector from C array:
  SEXP y_r;
  PROTECT(y_r = Rf_allocVector(REALSXP, static_cast<int>(size())));
  std::copy(y, y + size(), REAL(y_r));

  // 2. Slot for output
  SEXP dydt_r;

  // 3. Compute derivatives
  PROTECT(dydt_r = target(t, y_r));

  // 4. Copy derivatives from R vector into C array
  std::copy(REAL(dydt_r), REAL(dydt_r) + size(), dydt);
  
  // 5. Cleanup
  UNPROTECT(2);
}

