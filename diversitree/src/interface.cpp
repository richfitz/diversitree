#include <Rcpp.h>

#include "GslOdeBase.h"

#include "GslOdeR.h"
#include "GslOdeCompiled.h"

// Need to do this.
RCPP_EXPOSED_CLASS(GslOdeBase)

RCPP_MODULE(GslOde) {
  using namespace Rcpp;

  class_<GslOdeBase>("GslOdeBase")
    .property("size",      &GslOdeBase::size)
    .method("derivs",      &GslOdeBase::r_derivs)
    .method("run",         &GslOdeBase::r_run)
    .method("set_control", &GslOdeBase::r_set_control)
    ;

  class_<GslOdeR>("GslOdeR")
    .constructor<SEXP, SEXP, int>()
    .derives<GslOdeBase>("GslOdeBase")
    ;

  class_<GslOdeCompiled>("GslOdeCompiled")
    .constructor<SEXP, int>()
    .derives<GslOdeBase>("GslOdeBase")
    ;
}
