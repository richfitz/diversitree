#include <R.h>
#include <Rinternals.h>

#include "util-splines.h"

#include "hdr.h"

static double fn(double p, void *data) {
  HDRStruct dat = (HDRStruct) data;
  dt_spline *spline = dat->spline;
  return 
    dt_spline_eval1(spline, p + dat->w) - 
    dt_spline_eval1(spline, p);
}

SEXP hdr(SEXP x, SEXP y, SEXP alpha) {
  double xl, tol = 1e-8, a = REAL(alpha)[0];
  HDRStruct dat = (HDRStruct) R_alloc(1, sizeof(hdr_struct));
  SEXP ret;
  dat->spline = make_dt_spline(LENGTH(x), REAL(x), REAL(y), 0);
  dat->w      = 1-a;

  xl = RSRC_Brent_fmin(0, a, fn, (void *)dat, tol);
  PROTECT(ret = allocVector(REALSXP, 2));
  REAL(ret)[0] = dt_spline_eval1(dat->spline, xl);
  REAL(ret)[1] = dt_spline_eval1(dat->spline, xl+dat->w);
  UNPROTECT(1);

  return ret;
}

/* As it is non API, too, I've copied src/appl/fmin here.  Original
   copyright notice follows, and I claim no authorship over this
   code.  If this becomes API, I'll drop this */

/*
 *  R : A Computer Language for Statistical Data Analysis

 *  Copyright (C) 1997-2001  The R Development Core Team
 *  Copyright (C) 2003-2004  The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

double RSRC_Brent_fmin(double ax, double bx, double (*f)(double, void *),
		       void *info, double tol)
{
    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

/*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;/* -Wall */
    e = 0.;
    fx = (*f)(x, info);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

/*  main loop starts here ----------------------------------- */

    for(;;) {
	xm = (a + b) * .5;
	tol1 = eps * fabs(x) + tol3;
	t2 = tol1 * 2.;

	/* check stopping criterion */

	if (fabs(x - xm) <= t2 - (b - a) * .5) break;
	p = 0.;
	q = 0.;
	r = 0.;
	if (fabs(e) > tol1) { /* fit parabola */

	    r = (x - w) * (fx - fv);
	    q = (x - v) * (fx - fw);
	    p = (x - v) * q - (x - w) * r;
	    q = (q - r) * 2.;
	    if (q > 0.) p = -p; else q = -q;
	    r = e;
	    e = d;
	}

	if (fabs(p) >= fabs(q * .5 * r) ||
	    p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

	    if (x < xm) e = b - x; else e = a - x;
	    d = c * e;
	}
	else { /* a parabolic-interpolation step */

	    d = p / q;
	    u = x + d;

	    /* f must not be evaluated too close to ax or bx */

	    if (u - a < t2 || b - u < t2) {
		d = tol1;
		if (x >= xm) d = -d;
	    }
	}

	/* f must not be evaluated too close to x */

	if (fabs(d) >= tol1)
	    u = x + d;
	else if (d > 0.)
	    u = x + tol1;
	else
	    u = x - tol1;

	fu = (*f)(u, info);

	/*  update  a, b, v, w, and x */

	if (fu <= fx) {
	    if (u < x) b = x; else a = x;
	    v = w;    w = x;   x = u;
	    fv = fw; fw = fx; fx = fu;
	} else {
	    if (u < x) a = u; else b = u;
	    if (fu <= fw || w == x) {
		v = w; fv = fw;
		w = u; fw = fu;
	    } else if (fu <= fv || v == x || v == w) {
		v = u; fv = fu;
	    }
	}
    }
    /* end of main loop */

    return x;
}

