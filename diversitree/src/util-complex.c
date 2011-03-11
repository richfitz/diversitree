/* Utilities for working with complex numbers */
#include <R.h>
#include "util-complex.h"

/* Taken from ${R_SRC}/src/main/complex.c
   If we have C99 complex, then we can just use cexp(), which will
   probably be faster */
Rcomplex z_exp(Rcomplex z) {
  double expx;
  Rcomplex r;
  expx = exp(z.r);
  r.r = expx * cos(z.i);
  r.i = expx * sin(z.i);
  return r;
}

/* Addition */
Rcomplex z_plus(Rcomplex a, Rcomplex b) {
  Rcomplex c;
  c.r = a.r + b.r;
  c.i = a.i + b.i;
  return c;
}

/* Subtraction */
Rcomplex z_minus(Rcomplex a, Rcomplex b) {
  Rcomplex c;
  c.r = a.r - b.r;
  c.i = a.i - b.i;
  return c;
}

/* Multiplication */
Rcomplex z_times(Rcomplex a, Rcomplex b) {
  Rcomplex c;
  c.r = a.r * b.r - a.i * b.i;
  c.i = a.r * b.i + a.i * b.r;
  return c;
}

/* Division */
Rcomplex z_divide(Rcomplex a, Rcomplex b) {
  Rcomplex c;
  double ratio, den;
  double abr, abi;

  if( (abr = b.r) < 0)
    abr = - abr;
  if( (abi = b.i) < 0)
    abi = - abi;
  if( abr <= abi ) {
    ratio = b.r / b.i ;
    den = b.i * (1 + ratio*ratio);
    c.r = (a.r*ratio + a.i) / den;
    c.i = (a.i*ratio - a.r) / den;
  }
  else {
    ratio = b.i / b.r ;
    den = b.r * (1 + ratio*ratio);
    c.r = (a.r + a.i*ratio) / den;
    c.i = (a.i - a.r*ratio) / den;
  }

  return c;
}

/* Multiply by a scalar */
Rcomplex z_scal(Rcomplex a, double x) {
  Rcomplex c;
  c.r = a.r * x;
  c.i = a.i * x;
  return c;
}
