.onLoad <- function(libname, pkgname){
  .Call("set_sane_gsl_error_handling", PACKAGE="diversitree")
}
loadModule("diversitree", TRUE)
