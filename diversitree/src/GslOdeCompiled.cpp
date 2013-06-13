#include "GslOdeCompiled.h"
#include <R.h>

GslOdeCompiled::GslOdeCompiled(SEXP extPtr, size_t size) :
  GslOdeBase(size) {
  // This generates a warning that we may live with according to BDR:
  // https://stat.ethz.ch/pipermail/r-devel/2004-September/030792.html
  derivs_f = (deriv_func *) R_ExternalPtrAddr(extPtr);
}

void GslOdeCompiled::set_pars(SEXP pars_) {
  // Implicit short lifespan, as we don't copy (probably we
  // should...), but these are passed as read only elsewhere.
  // Declaring parameters as 
  //   double const * pars
  // should give us a mutable pointer to a constant vector (so that
  // pars can be reassigned, but the contents cannot be).  This
  // *should* work well here.
  pars = REAL(pars_);
}

void GslOdeCompiled::clear_pars() {
  pars = NULL;
}

void GslOdeCompiled::derivs(double t, const double y[], 
			    double dydt[]) {
  derivs_f(static_cast<int>(size()), t, pars, y, dydt);
}

