#include "GslOdeTime.h"
#include <R.h>
#include "util.h"

GslOdeTime::GslOdeTime(SEXP extPtr, size_t size, TimeMachine tm_) :
  GslOdeBase(size), tm(tm_) {
  derivs_f = get_deriv_func_from_ptr(extPtr);
}

void GslOdeTime::set_pars(SEXP pars_) {
  tm.set(Rcpp::as<std::vector<double> >(pars_));
}

void GslOdeTime::clear_pars() {
}

void GslOdeTime::derivs(double t, const double y[], 
			double dydt[]) {
  // TODO: An iterator approach might be better than the extra copy
  // that this involves?
  std::vector<double> pars = tm.get(t);
  derivs_f(static_cast<int>(size()), t, &pars[0], y, dydt);
}

