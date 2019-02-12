/* zeroin.f -- translated by f2c (version 20100827).
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

doublereal zeroin_(real *ax, real *bx, E_fp f, real *tol)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double r_sign(real *, real *);

    /* Local variables */
    static real a, b, c__, d__, e, p, q, r__, s, fa, fb, fc, xm, eps, tol1;


/*      a zero of the function  f(x)  is computed in the interval ax,bx . */

/*  input.. */

/*  ax     left endpoint of initial interval */
/*  bx     right endpoint of initial interval */
/*  f      function subprogram which evaluates f(x) for any x in */
/*         the interval  ax,bx */
/*  tol    desired length of the interval of uncertainty of the */
/*         final result ( .ge. 0.0) */


/*  output.. */

/*  zeroin abcissa approximating a zero of  f  in the interval ax,bx */


/*      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs */
/*  without  a  check.  zeroin  returns a zero  x  in the given interval */
/*  ax,bx  to within a tolerance  4*macheps*abs(x) + tol, where macheps */
/*  is the relative machine precision. */
/*      this function subprogram is a slightly  modified  translation  of */
/*  the algol 60 procedure  zero  given in  richard brent, algorithms for */
/*  minimization without derivatives, prentice - hall, inc. (1973). */



/*  compute eps, the relative machine precision */

    eps = 1.f;
L10:
    eps /= 2.f;
    tol1 = eps + 1.f;
    if (tol1 > 1.f) {
	goto L10;
    }

/* initialization */

    a = *ax;
    b = *bx;
    fa = (*f)(&a);
    fb = (*f)(&b);

/* begin step */

L20:
    c__ = a;
    fc = fa;
    d__ = b - a;
    e = d__;
L30:
    if (dabs(fc) >= dabs(fb)) {
	goto L40;
    }
    a = b;
    b = c__;
    c__ = a;
    fa = fb;
    fb = fc;
    fc = fa;

/* convergence test */

L40:
    tol1 = eps * 2.f * dabs(b) + *tol * .5f;
    xm = (c__ - b) * .5f;
    if (dabs(xm) <= tol1) {
	goto L90;
    }
    if (fb == 0.f) {
	goto L90;
    }

/* is bisection necessary */

    if (dabs(e) < tol1) {
	goto L70;
    }
    if (dabs(fa) <= dabs(fb)) {
	goto L70;
    }

/* is quadratic interpolation possible */

    if (a != c__) {
	goto L50;
    }

/* linear interpolation */

    s = fb / fa;
    p = xm * 2.f * s;
    q = 1.f - s;
    goto L60;

/* inverse quadratic interpolation */

L50:
    q = fa / fc;
    r__ = fb / fc;
    s = fb / fa;
    p = s * (xm * 2.f * q * (q - r__) - (b - a) * (r__ - 1.f));
    q = (q - 1.f) * (r__ - 1.f) * (s - 1.f);

/* adjust signs */

L60:
    if (p > 0.f) {
	q = -q;
    }
    p = dabs(p);

/* is interpolation acceptable */

    if (p * 2.f >= xm * 3.f * q - (r__1 = tol1 * q, dabs(r__1))) {
	goto L70;
    }
    if (p >= (r__1 = e * .5f * q, dabs(r__1))) {
	goto L70;
    }
    e = d__;
    d__ = p / q;
    goto L80;

/* bisection */

L70:
    d__ = xm;
    e = d__;

/* complete step */

L80:
    a = b;
    fa = fb;
    if (dabs(d__) > tol1) {
	b += d__;
    }
    if (dabs(d__) <= tol1) {
	b += r_sign(&tol1, &xm);
    }
    fb = (*f)(&b);
    if (fb * (fc / dabs(fc)) > 0.f) {
	goto L20;
    }
    goto L30;

/* done */

L90:
    ret_val = b;
    return ret_val;
} /* zeroin_ */

