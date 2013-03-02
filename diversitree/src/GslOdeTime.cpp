#include "GslOdeTime.h"
#include <R.h>

GslOdeTime::GslOdeTime(SEXP extPtr, int size, TimeMachine tm_) :
  GslOdeBase(size), tm(tm_) {
  // This generates a warning that we may live with according to BDR:
  // https://stat.ethz.ch/pipermail/r-devel/2004-September/030792.html
  derivs_f = (deriv_func *) R_ExternalPtrAddr(extPtr);
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
  derivs_f(size(), t, &pars[0], y, dydt);
}

