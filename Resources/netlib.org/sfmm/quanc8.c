/* quanc8.f -- translated by f2c (version 20100827).
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

/* Table of constant values */

static integer c__2 = 2;

/* Subroutine */ int quanc8_(E_fp fun, real *a, real *b, real *abserr, real *
	relerr, real *result, real *errest, integer *nofun, real *flag__)
{
    /* System generated locals */
    integer i__1;
    real r__1, r__2;

    /* Builtin functions */
    integer pow_ii(integer *, integer *);

    /* Local variables */
    static real f[16];
    static integer i__, j;
    static real x[16], f0, w0, w1, w2, w3, w4, x0;
    static integer nim, lev;
    static real area, cor11, temp, step, qnow, qdiff, fsave[240]	/* 
	    was [8][30] */;
    static integer nofin;
    static real qleft;
    static integer nomax;
    static real xsave[240]	/* was [8][30] */, stone, qprev;
    static integer levmin, levmax;
    static real qright[31], esterr, tolerr;
    static integer levout;



/*   estimate the integral of fun(x) from a to b */
/*   to a user provided tolerance. */
/*   an automatic adaptive routine based on */
/*   the 8-panel newton-cotes rule. */

/*   input .. */

/*   fun     the name of the integrand function subprogram fun(x). */
/*   a       the lower limit of integration. */
/*   b       the upper limit of integration.(b may be less than a.) */
/*   relerr  a relative error tolerance. (should be non-negative) */
/*   abserr  an absolute error tolerance. (should be non-negative) */

/*   output .. */

/*   result  an approximation to the integral hopefully satisfying the */
/*           least stringent of the two error tolerances. */
/*   errest  an estimate of the magnitude of the actual error. */
/*   nofun   the number of function values used in calculation of result. */
/*   flag    a reliability indicator.  if flag is zero, then result */
/*           probably satisfies the error tolerance.  if flag is */
/*           xxx.yyy , then  xxx = the number of intervals which have */
/*           not converged and  0.yyy = the fraction of the interval */
/*           left to do when the limit on  nofun  was approached. */


/*   ***   stage 1 ***   general initialization */
/*   set constants. */

    levmin = 1;
    levmax = 30;
    levout = 6;
    nomax = 5000;
    i__1 = levout + 1;
    nofin = nomax - (levmax - levout + pow_ii(&c__2, &i__1) << 3);

/*   trouble when nofun reaches nofin */

    w0 = .27908289241622575f;
    w1 = 1.6615167548500882f;
    w2 = -.26186948853615521f;
    w3 = 2.9618342151675483f;
    w4 = -1.2811287477954145f;

/*   initialize running sums to zero. */

    *flag__ = 0.f;
    *result = 0.f;
    cor11 = 0.f;
    *errest = 0.f;
    area = 0.f;
    *nofun = 0;
    if (*a == *b) {
	return 0;
    }

/*   ***   stage 2 ***   initialization for first interval */

    lev = 0;
    nim = 1;
    x0 = *a;
    x[15] = *b;
    qprev = 0.f;
    f0 = (*fun)(&x0);
    stone = (*b - *a) / 16.f;
    x[7] = (x0 + x[15]) / 2.f;
    x[3] = (x0 + x[7]) / 2.f;
    x[11] = (x[7] + x[15]) / 2.f;
    x[1] = (x0 + x[3]) / 2.f;
    x[5] = (x[3] + x[7]) / 2.f;
    x[9] = (x[7] + x[11]) / 2.f;
    x[13] = (x[11] + x[15]) / 2.f;
    for (j = 2; j <= 16; j += 2) {
	f[j - 1] = (*fun)(&x[j - 1]);
/* L25: */
    }
    *nofun = 9;

/*   ***   stage 3 ***   central calculation */
/*   requires qprev,x0,x2,x4,...,x16,f0,f2,f4,...,f16. */
/*   calculates x1,x3,...x15, f1,f3,...f15,qleft,qright,qnow,qdiff,area. */

L30:
    x[0] = (x0 + x[1]) / 2.f;
    f[0] = (*fun)(x);
    for (j = 3; j <= 15; j += 2) {
	x[j - 1] = (x[j - 2] + x[j]) / 2.f;
	f[j - 1] = (*fun)(&x[j - 1]);
/* L35: */
    }
    *nofun += 8;
    step = (x[15] - x0) / 16.f;
    qleft = (w0 * (f0 + f[7]) + w1 * (f[0] + f[6]) + w2 * (f[1] + f[5]) + w3 *
	     (f[2] + f[4]) + w4 * f[3]) * step;
    qright[lev] = (w0 * (f[7] + f[15]) + w1 * (f[8] + f[14]) + w2 * (f[9] + f[
	    13]) + w3 * (f[10] + f[12]) + w4 * f[11]) * step;
    qnow = qleft + qright[lev];
    qdiff = qnow - qprev;
    area += qdiff;

/*   ***   stage 4 *** interval convergence test */

    esterr = dabs(qdiff) / 1023.f;
/* Computing MAX */
    r__1 = *abserr, r__2 = *relerr * dabs(area);
    tolerr = dmax(r__1,r__2) * (step / stone);
    if (lev < levmin) {
	goto L50;
    }
    if (lev >= levmax) {
	goto L62;
    }
    if (*nofun > nofin) {
	goto L60;
    }
    if (esterr <= tolerr) {
	goto L70;
    }

/*   ***   stage 5   ***   no convergence */
/*   locate next interval. */

L50:
    nim <<= 1;
    ++lev;

/*   store right hand elements for future use. */

    for (i__ = 1; i__ <= 8; ++i__) {
	fsave[i__ + (lev << 3) - 9] = f[i__ + 7];
	xsave[i__ + (lev << 3) - 9] = x[i__ + 7];
/* L52: */
    }

/*   assemble left hand elements for immediate use. */

    qprev = qleft;
    for (i__ = 1; i__ <= 8; ++i__) {
	j = -i__;
	f[(j << 1) + 17] = f[j + 8];
	x[(j << 1) + 17] = x[j + 8];
/* L55: */
    }
    goto L30;

/*   ***   stage 6   ***   trouble section */
/*   number of function values is about to exceed limit. */

L60:
    nofin <<= 1;
    levmax = levout;
    *flag__ += (*b - x0) / (*b - *a);
    goto L70;

/*   current level is levmax. */

L62:
    *flag__ += 1.f;

/*   ***   stage 7   ***   interval converged */
/*   add contributions into running sums. */

L70:
    *result += qnow;
    *errest += esterr;
    cor11 += qdiff / 1023.f;

/*   locate next interval. */

L72:
    if (nim == nim / 2 << 1) {
	goto L75;
    }
    nim /= 2;
    --lev;
    goto L72;
L75:
    ++nim;
    if (lev <= 0) {
	goto L80;
    }

/*   assemble elements required for the next interval. */

    qprev = qright[lev - 1];
    x0 = x[15];
    f0 = f[15];
    for (i__ = 1; i__ <= 8; ++i__) {
	f[(i__ << 1) - 1] = fsave[i__ + (lev << 3) - 9];
	x[(i__ << 1) - 1] = xsave[i__ + (lev << 3) - 9];
/* L78: */
    }
    goto L30;

/*   ***   stage 8   ***   finalize and return */

L80:
    *result += cor11;

/*   make sure errest not less than roundoff level. */

    if (*errest == 0.f) {
	return 0;
    }
L82:
    temp = dabs(*result) + *errest;
    if (temp != dabs(*result)) {
	return 0;
    }
    *errest *= 2.f;
    goto L82;
} /* quanc8_ */

