/* spline.f -- translated by f2c (version 20100827).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Subroutine */ int spline_(integer *n, real *x, real *y, real *b, real *c__,
	 real *d__)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Local variables */
    static integer i__;
    static real t;
    static integer ib, nm1;


/*  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed */
/*  for a cubic interpolating spline */

/*    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3 */

/*    for  x(i) .le. x .le. x(i+1) */

/*  input.. */

/*    n = the number of data points or knots (n.ge.2) */
/*    x = the abscissas of the knots in strictly increasing order */
/*    y = the ordinates of the knots */

/*  output.. */

/*    b, c, d  = arrays of spline coefficients as defined above. */

/*  using  p  to denote differentiation, */

/*    y(i) = s(x(i)) */
/*    b(i) = sp(x(i)) */
/*    c(i) = spp(x(i))/2 */
/*    d(i) = sppp(x(i))/6  (derivative from the right) */

/*  the accompanying function subprogram  seval  can be used */
/*  to evaluate the spline. */



    /* Parameter adjustments */
    --d__;
    --c__;
    --b;
    --y;
    --x;

    /* Function Body */
    nm1 = *n - 1;
    if (*n < 2) {
	return 0;
    }
    if (*n < 3) {
	goto L50;
    }

/*  set up tridiagonal system */

/*  b = diagonal, d = offdiagonal, c = right hand side. */

    d__[1] = x[2] - x[1];
    c__[2] = (y[2] - y[1]) / d__[1];
    i__1 = nm1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	d__[i__] = x[i__ + 1] - x[i__];
	b[i__] = (d__[i__ - 1] + d__[i__]) * 2.f;
	c__[i__ + 1] = (y[i__ + 1] - y[i__]) / d__[i__];
	c__[i__] = c__[i__ + 1] - c__[i__];
/* L10: */
    }

/*  end conditions.  third derivatives at  x(1)  and  x(n) */
/*  obtained from divided differences */

    b[1] = -d__[1];
    b[*n] = -d__[*n - 1];
    c__[1] = 0.f;
    c__[*n] = 0.f;
    if (*n == 3) {
	goto L15;
    }
    c__[1] = c__[3] / (x[4] - x[2]) - c__[2] / (x[3] - x[1]);
    c__[*n] = c__[*n - 1] / (x[*n] - x[*n - 2]) - c__[*n - 2] / (x[*n - 1] - 
	    x[*n - 3]);
/* Computing 2nd power */
    r__1 = d__[1];
    c__[1] = c__[1] * (r__1 * r__1) / (x[4] - x[1]);
/* Computing 2nd power */
    r__1 = d__[*n - 1];
    c__[*n] = -c__[*n] * (r__1 * r__1) / (x[*n] - x[*n - 3]);

/*  forward elimination */

L15:
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	t = d__[i__ - 1] / b[i__ - 1];
	b[i__] -= t * d__[i__ - 1];
	c__[i__] -= t * c__[i__ - 1];
/* L20: */
    }

/*  back substitution */

    c__[*n] /= b[*n];
    i__1 = nm1;
    for (ib = 1; ib <= i__1; ++ib) {
	i__ = *n - ib;
	c__[i__] = (c__[i__] - d__[i__] * c__[i__ + 1]) / b[i__];
/* L30: */
    }

/*  c(i) is now the sigma(i) of the text */

/*  compute polynomial coefficients */

    b[*n] = (y[*n] - y[nm1]) / d__[nm1] + d__[nm1] * (c__[nm1] + c__[*n] * 
	    2.f);
    i__1 = nm1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b[i__] = (y[i__ + 1] - y[i__]) / d__[i__] - d__[i__] * (c__[i__ + 1] 
		+ c__[i__] * 2.f);
	d__[i__] = (c__[i__ + 1] - c__[i__]) / d__[i__];
	c__[i__] *= 3.f;
/* L40: */
    }
    c__[*n] *= 3.f;
    d__[*n] = d__[*n - 1];
    return 0;

L50:
    b[1] = (y[2] - y[1]) / (x[2] - x[1]);
    c__[1] = 0.f;
    d__[1] = 0.f;
    b[2] = b[1];
    c__[2] = 0.f;
    d__[2] = 0.f;
    return 0;
} /* spline_ */

