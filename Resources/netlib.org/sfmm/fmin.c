/* fmin.f -- translated by f2c (version 20100827).
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

doublereal fmin_(real *ax, real *bx, E_fp f, real *tol)
{
    /* System generated locals */
    real ret_val, r__1;

    /* Builtin functions */
    double sqrt(doublereal), r_sign(real *, real *);

    /* Local variables */
    static real a, b, c__, d__, e, p, q, r__, u, v, w, x, fu, fv, fw, fx, xm, 
	    eps, tol1, tol2;


/*      an approximation  x  to the point where  f  attains a minimum  on */
/*  the interval  (ax,bx)  is determined. */


/*  input.. */

/*  ax    left endpoint of initial interval */
/*  bx    right endpoint of initial interval */
/*  f     function subprogram which evaluates  f(x)  for any  x */
/*        in the interval  (ax,bx) */
/*  tol   desired length of the interval of uncertainty of the final */
/*        result ( .ge. 0.0) */


/*  output.. */

/*  fmin  abcissa approximating the point where  f  attains a minimum */


/*      the method used is a combination of  golden  section  search  and */
/*  successive parabolic interpolation.  convergence is never much slower */
/*  than  that  for  a  fibonacci search.  if  f  has a continuous second */
/*  derivative which is positive at the minimum (which is not  at  ax  or */
/*  bx),  then  convergence  is  superlinear, and usually of the order of */
/*  about  1.324.... */
/*      the function  f  is never evaluated at two points closer together */
/*  than  eps*abs(fmin) + (tol/3), where eps is  approximately the square */
/*  root  of  the  relative  machine  precision.   if   f   is a unimodal */
/*  function and the computed values of   f   are  always  unimodal  when */
/*  separated by at least  eps*abs(x) + (tol/3), then  fmin  approximates */
/*  the abcissa of the global minimum of  f  on the interval  ax,bx  with */
/*  an error less than  3*eps*abs(fmin) + tol.  if   f   is not unimodal, */
/*  then fmin may approximate a local, but perhaps non-global, minimum to */
/*  the same accuracy. */
/*      this function subprogram is a slightly modified  version  of  the */
/*  algol  60 procedure  localmin  given in richard brent, algorithms for */
/*  minimization without derivatives, prentice - hall, inc. (1973). */



/*  c is the squared inverse of the golden ratio */

    c__ = (3.f - sqrt(5.f)) * .5f;

/*  eps is approximately the square root of the relative machine */
/*  precision. */

    eps = 1.f;
L10:
    eps /= 2.f;
    tol1 = eps + 1.f;
    if (tol1 > 1.f) {
	goto L10;
    }
    eps = sqrt(eps);

/*  initialization */

    a = *ax;
    b = *bx;
    v = a + c__ * (b - a);
    w = v;
    x = v;
    e = 0.f;
    fx = (*f)(&x);
    fv = fx;
    fw = fx;

/*  main loop starts here */

L20:
    xm = (a + b) * .5f;
    tol1 = eps * dabs(x) + *tol / 3.f;
    tol2 = tol1 * 2.f;

/*  check stopping criterion */

    if ((r__1 = x - xm, dabs(r__1)) <= tol2 - (b - a) * .5f) {
	goto L90;
    }

/* is golden-section necessary */

    if (dabs(e) <= tol1) {
	goto L40;
    }

/*  fit parabola */

    r__ = (x - w) * (fx - fv);
    q = (x - v) * (fx - fw);
    p = (x - v) * q - (x - w) * r__;
    q = (q - r__) * 2.f;
    if (q > 0.f) {
	p = -p;
    }
    q = dabs(q);
    r__ = e;
    e = d__;

/*  is parabola acceptable */

/* L30: */
    if (dabs(p) >= (r__1 = q * .5f * r__, dabs(r__1))) {
	goto L40;
    }
    if (p <= q * (a - x)) {
	goto L40;
    }
    if (p >= q * (b - x)) {
	goto L40;
    }

/*  a parabolic interpolation step */

    d__ = p / q;
    u = x + d__;

/*  f must not be evaluated too close to ax or bx */

    if (u - a < tol2) {
	r__1 = xm - x;
	d__ = r_sign(&tol1, &r__1);
    }
    if (b - u < tol2) {
	r__1 = xm - x;
	d__ = r_sign(&tol1, &r__1);
    }
    goto L50;

/*  a golden-section step */

L40:
    if (x >= xm) {
	e = a - x;
    }
    if (x < xm) {
	e = b - x;
    }
    d__ = c__ * e;

/*  f must not be evaluated too close to x */

L50:
    if (dabs(d__) >= tol1) {
	u = x + d__;
    }
    if (dabs(d__) < tol1) {
	u = x + r_sign(&tol1, &d__);
    }
    fu = (*f)(&u);

/*  update  a, b, v, w, and x */

    if (fu > fx) {
	goto L60;
    }
    if (u >= x) {
	a = x;
    }
    if (u < x) {
	b = x;
    }
    v = w;
    fv = fw;
    w = x;
    fw = fx;
    x = u;
    fx = fu;
    goto L20;
L60:
    if (u < x) {
	a = u;
    }
    if (u >= x) {
	b = u;
    }
    if (fu <= fw) {
	goto L70;
    }
    if (w == x) {
	goto L70;
    }
    if (fu <= fv) {
	goto L80;
    }
    if (v == x) {
	goto L80;
    }
    if (v == w) {
	goto L80;
    }
    goto L20;
L70:
    v = w;
    fv = fw;
    w = u;
    fw = fu;
    goto L20;
L80:
    v = u;
    fv = fu;
    goto L20;

/*  end of main loop */

L90:
    ret_val = x;
    return ret_val;
} /* fmin_ */

