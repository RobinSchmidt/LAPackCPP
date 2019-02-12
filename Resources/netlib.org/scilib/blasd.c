/* blasd.f -- translated by f2c (version 20100827).
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

static doublereal c_b57 = 1.;
static doublecomplex c_b109 = {1.,0.};
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__9 = 9;

/* -------------------------------------------------------------     ************ */
/*                                                                     SASUM */
/*                                                                  ************ */
doublereal sasum_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4, d__5, d__6;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, m, ix, mp1;
    static doublereal dtemp;


/*     TAKES THE SUM OF THE ABSOLUTE VALUES. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += (d__1 = dx[ix], abs(d__1));
	ix += *incx;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        CODE FOR INCREMENT EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += (d__1 = dx[i__], abs(d__1));
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 6) {
	dtemp = dtemp + (d__1 = dx[i__], abs(d__1)) + (d__2 = dx[i__ + 1], 
		abs(d__2)) + (d__3 = dx[i__ + 2], abs(d__3)) + (d__4 = dx[i__ 
		+ 3], abs(d__4)) + (d__5 = dx[i__ + 4], abs(d__5)) + (d__6 = 
		dx[i__ + 5], abs(d__6));
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* sasum_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SAXPY */
/*                                                                  ************ */
/* Subroutine */ int saxpy_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx, doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, m, ix, iy, mp1;


/*     CONSTANT TIMES A VECTOR PLUS A VECTOR. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] += *da * dx[i__];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	dy[i__] += *da * dx[i__];
	dy[i__ + 1] += *da * dx[i__ + 1];
	dy[i__ + 2] += *da * dx[i__ + 2];
	dy[i__ + 3] += *da * dx[i__ + 3];
/* L50: */
    }
    return 0;
} /* saxpy_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SCOPY */
/*                                                                  ************ */
/* Subroutine */ int scopy_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, m, ix, iy, mp1;


/*     COPIES A VECTOR, X, TO A VECTOR, Y. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] = dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] = dx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	dy[i__] = dx[i__];
	dy[i__ + 1] = dx[i__ + 1];
	dy[i__ + 2] = dx[i__ + 2];
	dy[i__ + 3] = dx[i__ + 3];
	dy[i__ + 4] = dx[i__ + 4];
	dy[i__ + 5] = dx[i__ + 5];
	dy[i__ + 6] = dx[i__ + 6];
/* L50: */
    }
    return 0;
} /* scopy_ */

/* -------------------------------------------------------------     ************ */
/*                                                                      SDOT */
/*                                                                  ************ */
doublereal sdot_(integer *n, doublereal *dx, integer *incx, doublereal *dy, 
	integer *incy)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static doublereal dtemp;


/*     FORMS THE DOT PRODUCT OF TWO VECTORS. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[i__] * dy[i__];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
		i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 
		4] * dy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* sdot_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SMACH */
/*                                                                  ************ */
doublereal smach_(integer *job)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal s, eps, huge__, tiny;


/*     SMACH COMPUTES MACHINE PARAMETERS OF FLOATING POINT */
/*     ARITHMETIC FOR USE IN TESTING ONLY.  NOT REQUIRED BY */
/*     LINPACK PROPER. */

/*     IF TROUBLE WITH AUTOMATIC COMPUTATION OF THESE QUANTITIES, */
/*     THEY CAN BE SET BY DIRECT ASSIGNMENT STATEMENTS. */
/*     ASSUME THE COMPUTER HAS */

/*        B = BASE OF ARITHMETIC */
/*        T = NUMBER OF BASE  B  DIGITS */
/*        L = SMALLEST POSSIBLE EXPONENT */
/*        U = LARGEST POSSIBLE EXPONENT */

/*     THEN */

/*        EPS = B**(1-T) */
/*        TINY = 100.0*B**(-L+T) */
/*        HUGE = 0.01*B**(U-T) */

/*     DMACH SAME AS SMACH EXCEPT T, L, U APPLY TO */
/*     DOUBLE PRECISION. */

/*     CMACH SAME AS SMACH EXCEPT IF COMPLEX DIVISION */
/*     IS DONE BY */

/*        1/(X+I*Y) = (X-I*Y)/(X**2+Y**2) */

/*     THEN */

/*        TINY = SQRT(TINY) */
/*        HUGE = SQRT(HUGE) */


/*     JOB IS 1, 2 OR 3 FOR EPSILON, TINY AND HUGE, RESPECTIVELY. */


    eps = 1.;
L10:
    eps /= 2.;
    s = eps + 1.;
    if (s > 1.) {
	goto L10;
    }
    eps *= 2.;

    s = 1.;
L20:
    tiny = s;
    s /= 16.;
    if (s * 1.f != 0.) {
	goto L20;
    }
    tiny = tiny / eps * 100.f;
    huge__ = 1. / tiny;

    if (*job == 1) {
	ret_val = eps;
    }
    if (*job == 2) {
	ret_val = tiny;
    }
    if (*job == 3) {
	ret_val = huge__;
    }
    return ret_val;
} /* smach_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SNRM2 */
/*                                                                  ************ */
doublereal snrm2_(integer *n, doublereal *dx, integer *incx)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal cutlo = 8.232e-11;
    static doublereal cuthi = 1.304e19;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, nn;
    static doublereal sum, xmax;
    static integer next;
    static doublereal hitest;

    /* Assigned format variables */
    static char *next_fmt;

    /* Parameter adjustments */
    --dx;

    /* Function Body */

/*     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE */
/*     INCREMENT INCX . */
/*     IF    N .LE. 0 RETURN WITH RESULT = 0. */
/*     IF N .GE. 1 THEN INCX MUST BE .GE. 1 */

/*           C.L.LAWSON, 1978 JAN 08 */

/*     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE */
/*     HOPEFULLY APPLICABLE TO ALL MACHINES. */
/*         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES. */
/*         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES. */
/*     WHERE */
/*         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1. */
/*         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT) */
/*         V   = LARGEST  NO.            (OVERFLOW  LIMIT) */

/*     BRIEF OUTLINE OF ALGORITHM.. */

/*     PHASE 1    SCANS ZERO COMPONENTS. */
/*     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO */
/*     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO */
/*     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M */
/*     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX. */

/*     VALUES FOR CUTLO AND CUTHI.. */
/*     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER */
/*     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS.. */
/*     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE */
/*                   UNIVAC AND DEC AT 2**(-103) */
/*                   THUS CUTLO = 2**(-51) = 4.44089E-16 */
/*     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC. */
/*                   THUS CUTHI = 2**(63.5) = 1.30438E19 */
/*     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC. */
/*                   THUS CUTLO = 2**(-33.5) = 8.23181D-11 */
/*     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19 */
/*     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 / */
/*     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 / */

    if (*n > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    nn = *n * *incx;
/*                    BEGIN MAIN LOOP */
    i__ = 1;
L20:
    switch (next) {
	case 0: goto L30;
	case 1: goto L50;
	case 2: goto L70;
	case 3: goto L110;
    }
L30:
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L85;
    }
    next = 1;
    next_fmt = fmt_50;
    xmax = zero;

/*                        PHASE 1.  SUM IS ZERO */

L50:
    if (dx[i__] == zero) {
	goto L200;
    }
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L85;
    }

/*                                PREPARE FOR PHASE 2. */
    next = 2;
    next_fmt = fmt_70;
    goto L105;

/*                                PREPARE FOR PHASE 4. */

L100:
    i__ = j;
    next = 3;
    next_fmt = fmt_110;
    sum = sum / dx[i__] / dx[i__];
L105:
    xmax = (d__1 = dx[i__], abs(d__1));
    goto L115;

/*                   PHASE 2.  SUM IS SMALL. */
/*                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW. */

L70:
    if ((d__1 = dx[i__], abs(d__1)) > cutlo) {
	goto L75;
    }

/*                     COMMON CODE FOR PHASES 2 AND 4. */
/*                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. */

L110:
    if ((d__1 = dx[i__], abs(d__1)) <= xmax) {
	goto L115;
    }
/* Computing 2nd power */
    d__1 = xmax / dx[i__];
    sum = one + sum * (d__1 * d__1);
    xmax = (d__1 = dx[i__], abs(d__1));
    goto L200;

L115:
/* Computing 2nd power */
    d__1 = dx[i__] / xmax;
    sum += d__1 * d__1;
    goto L200;


/*                  PREPARE FOR PHASE 3. */

L75:
    sum = sum * xmax * xmax;


/*     FOR REAL OR D.P. SET HITEST = CUTHI/N */
/*     FOR COMPLEX      SET HITEST = CUTHI/(2*N) */

L85:
    hitest = cuthi / (real) (*n);

/*                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING. */

    i__1 = nn;
    i__2 = *incx;
    for (j = i__; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
	if ((d__1 = dx[j], abs(d__1)) >= hitest) {
	    goto L100;
	}
/* L95: */
/* Computing 2nd power */
	d__1 = dx[j];
	sum += d__1 * d__1;
    }
    ret_val = sqrt(sum);
    goto L300;

L200:
    i__ += *incx;
    if (i__ <= nn) {
	goto L20;
    }

/*              END OF MAIN LOOP. */

/*              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING. */

    ret_val = xmax * sqrt(sum);
L300:
    return ret_val;
} /* snrm2_ */

/* -------------------------------------------------------------     ************ */
/*                                                                      SROT */
/*                                                                  ************ */
/* Subroutine */ int srot_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, ix, iy;
    static doublereal dtemp;


/*     APPLIES A PLANE ROTATION. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL */
/*         TO 1 */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = *c__ * dx[ix] + *s * dy[iy];
	dy[iy] = *c__ * dy[iy] - *s * dx[ix];
	dx[ix] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       CODE FOR BOTH INCREMENTS EQUAL TO 1 */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = *c__ * dx[i__] + *s * dy[i__];
	dy[i__] = *c__ * dy[i__] - *s * dx[i__];
	dx[i__] = dtemp;
/* L30: */
    }
    return 0;
} /* srot_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SROTG */
/*                                                                  ************ */
/* Subroutine */ int srotg_(doublereal *da, doublereal *db, doublereal *c__, 
	doublereal *s)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal r__, z__, roe, scale;


/*     CONSTRUCT GIVENS PLANE ROTATION. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    roe = *db;
    if (abs(*da) > abs(*db)) {
	roe = *da;
    }
    scale = abs(*da) + abs(*db);
    if (scale != 0.) {
	goto L10;
    }
    *c__ = 1.;
    *s = 0.;
    r__ = 0.;
    goto L20;
L10:
/* Computing 2nd power */
    d__1 = *da / scale;
/* Computing 2nd power */
    d__2 = *db / scale;
    r__ = scale * sqrt(d__1 * d__1 + d__2 * d__2);
    r__ = d_sign(&c_b57, &roe) * r__;
    *c__ = *da / r__;
    *s = *db / r__;
L20:
    z__ = 1.;
    if (abs(*da) > abs(*db)) {
	z__ = *s;
    }
    if (abs(*db) >= abs(*da) && *c__ != 0.) {
	z__ = 1. / *c__;
    }
    *da = r__;
    *db = z__;
    return 0;
} /* srotg_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SSCAL */
/*                                                                  ************ */
/* Subroutine */ int sscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, m, ix, mp1;


/*     SCALES A VECTOR BY A CONSTANT. */
/*     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dx[ix] = *da * dx[ix];
	ix += *incx;
/* L10: */
    }
    return 0;

/*        CODE FOR INCREMENT EQUAL TO 1 */


/*        CLEAN-UP LOOP */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dx[i__] = *da * dx[i__];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	dx[i__] = *da * dx[i__];
	dx[i__ + 1] = *da * dx[i__ + 1];
	dx[i__ + 2] = *da * dx[i__ + 2];
	dx[i__ + 3] = *da * dx[i__ + 3];
	dx[i__ + 4] = *da * dx[i__ + 4];
/* L50: */
    }
    return 0;
} /* sscal_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SSWAP */
/*                                                                  ************ */
/* Subroutine */ int sswap_(integer *n, doublereal *dx, integer *incx, 
	doublereal *dy, integer *incy)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, m, ix, iy, mp1;
    static doublereal dtemp;


/*     INTERCHANGES TWO VECTORS. */
/*     USES UNROLLED LOOPS FOR INCREMENTS EQUAL ONE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL */
/*         TO 1 */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = dx[ix];
	dx[ix] = dy[iy];
	dy[iy] = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       CODE FOR BOTH INCREMENTS EQUAL TO 1 */


/*       CLEAN-UP LOOP */

L20:
    m = *n % 3;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp = dx[i__];
	dx[i__] = dy[i__];
	dy[i__] = dtemp;
/* L30: */
    }
    if (*n < 3) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 3) {
	dtemp = dx[i__];
	dx[i__] = dy[i__];
	dy[i__] = dtemp;
	dtemp = dx[i__ + 1];
	dx[i__ + 1] = dy[i__ + 1];
	dy[i__ + 1] = dtemp;
	dtemp = dx[i__ + 2];
	dx[i__ + 2] = dy[i__ + 2];
	dy[i__ + 2] = dtemp;
/* L50: */
    }
    return 0;
} /* sswap_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     ISAMAX */
/*                                                                  ************ */
integer isamax_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, ix;
    static doublereal dmax__;


/*     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    dmax__ = (d__1 = dx[ix], abs(d__1));
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if ((d__1 = dx[ix], abs(d__1)) <= dmax__) {
	    goto L5;
	}
	ret_val = i__;
	dmax__ = (d__1 = dx[ix], abs(d__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENT EQUAL TO 1 */

L20:
    dmax__ = abs(dx[1]);
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if ((d__1 = dx[i__], abs(d__1)) <= dmax__) {
	    goto L30;
	}
	ret_val = i__;
	dmax__ = (d__1 = dx[i__], abs(d__1));
L30:
	;
    }
    return ret_val;
} /* isamax_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SPDOT */
/*                                                                  ************ */
doublereal spdot_(integer *n, doublereal *sy, integer *index, doublereal *sx)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__;
    static doublereal sum;


    /* Parameter adjustments */
    --sx;
    --index;
    --sy;

    /* Function Body */
    sum = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	sum += sy[index[i__]] * sx[i__];
    }
    ret_val = sum;
    return ret_val;
} /* spdot_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SPAXPY */
/*                                                                  ************ */
/* Subroutine */ int spaxpy_(integer *n, doublereal *sa, doublereal *sx, 
	doublereal *sy, integer *index)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;


    /* Parameter adjustments */
    --index;
    --sy;
    --sx;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* L60: */
	sy[index[j]] = *sa * sx[j] + sy[index[j]];
    }
    return 0;
} /* spaxpy_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     CAXPY */
/*                                                                  ************ */
/* Subroutine */ int caxpy_(integer *n, doublecomplex *ca, doublecomplex *cx, 
	integer *incx, doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, ix, iy;


/*     CONSTANT TIMES A VECTOR PLUS A VECTOR. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if ((d__1 = ca->r, abs(d__1)) + (d__2 = d_imag(ca), abs(d__2)) == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iy;
	i__3 = iy;
	i__4 = ix;
	z__2.r = ca->r * cx[i__4].r - ca->i * cx[i__4].i, z__2.i = ca->r * cx[
		i__4].i + ca->i * cx[i__4].r;
	z__1.r = cy[i__3].r + z__2.r, z__1.i = cy[i__3].i + z__2.i;
	cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	i__4 = i__;
	z__2.r = ca->r * cx[i__4].r - ca->i * cx[i__4].i, z__2.i = ca->r * cx[
		i__4].i + ca->i * cx[i__4].r;
	z__1.r = cy[i__3].r + z__2.r, z__1.i = cy[i__3].i + z__2.i;
	cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
/* L30: */
    }
    return 0;
} /* caxpy_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     CCOPY */
/*                                                                  ************ */
/* Subroutine */ int ccopy_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, ix, iy;


/*     COPIES A VECTOR, X, TO A VECTOR, Y. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = iy;
	i__3 = ix;
	cy[i__2].r = cx[i__3].r, cy[i__2].i = cx[i__3].i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	cy[i__2].r = cx[i__3].r, cy[i__2].i = cx[i__3].i;
/* L30: */
    }
    return 0;
} /* ccopy_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     CDOTC */
/*                                                                  ************ */
/* Double Complex */ VOID cdotc_(doublecomplex * ret_val, integer *n, 
	doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ctemp;


/*     FORMS THE DOT PRODUCT OF TWO VECTORS, CONJUGATING THE FIRST */
/*     VECTOR. */
/*     JACK DONGARRA, LINPACK,  3/11/78. */


    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    ctemp.r = 0., ctemp.i = 0.;
     ret_val->r = 0.,  ret_val->i = 0.;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d_cnjg(&z__3, &cx[ix]);
	i__2 = iy;
	z__2.r = z__3.r * cy[i__2].r - z__3.i * cy[i__2].i, z__2.i = z__3.r * 
		cy[i__2].i + z__3.i * cy[i__2].r;
	z__1.r = ctemp.r + z__2.r, z__1.i = ctemp.i + z__2.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
    return ;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d_cnjg(&z__3, &cx[i__]);
	i__2 = i__;
	z__2.r = z__3.r * cy[i__2].r - z__3.i * cy[i__2].i, z__2.i = z__3.r * 
		cy[i__2].i + z__3.i * cy[i__2].r;
	z__1.r = ctemp.r + z__2.r, z__1.i = ctemp.i + z__2.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
/* L30: */
    }
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
    return ;
} /* cdotc_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     CDOTU */
/*                                                                  ************ */
/* Double Complex */ VOID cdotu_(doublecomplex * ret_val, integer *n, 
	doublecomplex *cx, integer *incx, doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ctemp;


/*     FORMS THE DOT PRODUCT OF TWO VECTORS. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    ctemp.r = 0., ctemp.i = 0.;
     ret_val->r = 0.,  ret_val->i = 0.;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*          NOT EQUAL TO 1 */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ix;
	i__3 = iy;
	z__2.r = cx[i__2].r * cy[i__3].r - cx[i__2].i * cy[i__3].i, z__2.i = 
		cx[i__2].r * cy[i__3].i + cx[i__2].i * cy[i__3].r;
	z__1.r = ctemp.r + z__2.r, z__1.i = ctemp.i + z__2.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
    return ;

/*        CODE FOR BOTH INCREMENTS EQUAL TO 1 */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	z__2.r = cx[i__2].r * cy[i__3].r - cx[i__2].i * cy[i__3].i, z__2.i = 
		cx[i__2].r * cy[i__3].i + cx[i__2].i * cy[i__3].r;
	z__1.r = ctemp.r + z__2.r, z__1.i = ctemp.i + z__2.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
/* L30: */
    }
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
    return ;
} /* cdotu_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     CMACH */
/*                                                                  ************ */
doublereal cmach_(integer *job)
{
    /* System generated locals */
    doublereal ret_val;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal s, eps, huge__, tiny;


/*     SMACH COMPUTES MACHINE PARAMETERS OF FLOATING POINT */
/*     ARITHMETIC FOR USE IN TESTING ONLY.  NOT REQUIRED BY */
/*     LINPACK PROPER. */

/*     IF TROUBLE WITH AUTOMATIC COMPUTATION OF THESE QUANTITIES, */
/*     THEY CAN BE SET BY DIRECT ASSIGNMENT STATEMENTS. */
/*     ASSUME THE COMPUTER HAS */

/*        B = BASE OF ARITHMETIC */
/*        T = NUMBER OF BASE  B  DIGITS */
/*        L = SMALLEST POSSIBLE EXPONENT */
/*        U = LARGEST POSSIBLE EXPONENT */

/*     THEN */

/*        EPS = B**(1-T) */
/*        TINY = 100.0*B**(-L+T) */
/*        HUGE = 0.01*B**(U-T) */

/*     DMACH SAME AS SMACH EXCEPT T, L, U APPLY TO */
/*     DOUBLE PRECISION. */

/*     CMACH SAME AS SMACH EXCEPT IF COMPLEX DIVISION */
/*     IS DONE BY */

/*        1/(X+I*Y) = (X-I*Y)/(X**2+Y**2) */

/*     THEN */

/*        TINY = SQRT(TINY) */
/*        HUGE = SQRT(HUGE) */


/*     JOB IS 1, 2 OR 3 FOR EPSILON, TINY AND HUGE, RESPECTIVELY. */



    eps = 1.;
L10:
    eps /= 2.;
    s = eps + 1.f;
    if (s > 1.) {
	goto L10;
    }
    eps *= 2.;
    ret_val = eps;
    if (*job == 1) {
	return ret_val;
    }

    s = 1.;
L20:
    tiny = s;
    s /= 16.;
    if (s * 1. != 0.) {
	goto L20;
    }
    tiny = tiny / eps * 100.;
    z__2.r = tiny, z__2.i = 0.;
    z_div(&z__1, &c_b109, &z__2);
    s = z__1.r;
    if (s != 1. / tiny) {
	tiny = sqrt(tiny);
    }
    huge__ = 1. / tiny;
    if (*job == 1) {
	ret_val = eps;
    }
    if (*job == 2) {
	ret_val = tiny;
    }
    if (*job == 3) {
	ret_val = huge__;
    }
    return ret_val;
} /* cmach_ */

/* -------------------------------------------------------------     ************ */
/*                                                                      CROT */
/*                                                                  ************ */
/* Subroutine */ int crot_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublecomplex *cc, doublecomplex *
	cs)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ctemp;


/*     APPLIES A PLANE ROTATION. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL */
/*         TO 1. */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ix;
	z__2.r = cc->r * cx[i__2].r - cc->i * cx[i__2].i, z__2.i = cc->r * cx[
		i__2].i + cc->i * cx[i__2].r;
	i__3 = iy;
	z__3.r = cs->r * cy[i__3].r - cs->i * cy[i__3].i, z__3.i = cs->r * cy[
		i__3].i + cs->i * cy[i__3].r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	i__2 = iy;
	i__3 = iy;
	z__2.r = cc->r * cy[i__3].r - cc->i * cy[i__3].i, z__2.i = cc->r * cy[
		i__3].i + cc->i * cy[i__3].r;
	d_cnjg(&z__4, cs);
	i__4 = ix;
	z__3.r = z__4.r * cx[i__4].r - z__4.i * cx[i__4].i, z__3.i = z__4.r * 
		cx[i__4].i + z__4.i * cx[i__4].r;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
	i__2 = ix;
	cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*       CODE FOR BOTH INCREMENTS EQUAL TO 1 */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	z__2.r = cc->r * cx[i__2].r - cc->i * cx[i__2].i, z__2.i = cc->r * cx[
		i__2].i + cc->i * cx[i__2].r;
	i__3 = i__;
	z__3.r = cs->r * cy[i__3].r - cs->i * cy[i__3].i, z__3.i = cs->r * cy[
		i__3].i + cs->i * cy[i__3].r;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	i__2 = i__;
	i__3 = i__;
	z__2.r = cc->r * cy[i__3].r - cc->i * cy[i__3].i, z__2.i = cc->r * cy[
		i__3].i + cc->i * cy[i__3].r;
	d_cnjg(&z__4, cs);
	i__4 = i__;
	z__3.r = z__4.r * cx[i__4].r - z__4.i * cx[i__4].i, z__3.i = z__4.r * 
		cx[i__4].i + z__4.i * cx[i__4].r;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
	i__2 = i__;
	cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
/* L30: */
    }
    return 0;
} /* crot_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     CROTG */
/*                                                                  ************ */
/* Subroutine */ int crotg_(doublecomplex *ca, doublecomplex *cb, doublereal *
	c__, doublecomplex *s)
{
    /* System generated locals */
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double z_abs(doublecomplex *), sqrt(doublereal);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal norm;
    static doublecomplex alpha;
    static doublereal scale;

    if (z_abs(ca) != 0.) {
	goto L10;
    }
    *c__ = 0.;
    s->r = 1., s->i = 0.;
    ca->r = cb->r, ca->i = cb->i;
    goto L20;
L10:
    scale = z_abs(ca) + z_abs(cb);
    z__1.r = ca->r / scale, z__1.i = ca->i / scale;
/* Computing 2nd power */
    d__1 = z_abs(&z__1);
    z__2.r = cb->r / scale, z__2.i = cb->i / scale;
/* Computing 2nd power */
    d__2 = z_abs(&z__2);
    norm = scale * sqrt(d__1 * d__1 + d__2 * d__2);
    d__1 = z_abs(ca);
    z__1.r = ca->r / d__1, z__1.i = ca->i / d__1;
    alpha.r = z__1.r, alpha.i = z__1.i;
    *c__ = z_abs(ca) / norm;
    d_cnjg(&z__3, cb);
    z__2.r = alpha.r * z__3.r - alpha.i * z__3.i, z__2.i = alpha.r * z__3.i + 
	    alpha.i * z__3.r;
    z__1.r = z__2.r / norm, z__1.i = z__2.i / norm;
    s->r = z__1.r, s->i = z__1.i;
    z__1.r = norm * alpha.r, z__1.i = norm * alpha.i;
    ca->r = z__1.r, ca->i = z__1.i;
L20:
    return 0;
} /* crotg_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     CSCAL */
/*                                                                  ************ */
/* Subroutine */ int cscal_(integer *n, doublecomplex *ca, doublecomplex *cx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, ix;


/*     SCALES A VECTOR BY A CONSTANT. */
/*     JACK DONGARRA, LINPACK,  3/11/78. */


    /* Parameter adjustments */
    --cx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1. */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ix;
	i__3 = ix;
	z__1.r = ca->r * cx[i__3].r - ca->i * cx[i__3].i, z__1.i = ca->r * cx[
		i__3].i + ca->i * cx[i__3].r;
	cx[i__2].r = z__1.r, cx[i__2].i = z__1.i;
	ix += *incx;
/* L10: */
    }
    return 0;

/*        CODE FOR INCREMENT EQUAL TO 1. */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	z__1.r = ca->r * cx[i__3].r - ca->i * cx[i__3].i, z__1.i = ca->r * cx[
		i__3].i + ca->i * cx[i__3].r;
	cx[i__2].r = z__1.r, cx[i__2].i = z__1.i;
/* L30: */
    }
    return 0;
} /* cscal_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     CSSCAL */
/*                                                                  ************ */
/* Subroutine */ int csscal_(integer *n, doublereal *sa, doublecomplex *cx, 
	integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, ix;


/*     SCALES A COMPLEX VECTOR BY A REAL CONSTANT. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --cx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1. */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ix;
	i__3 = ix;
	d__1 = *sa * cx[i__3].r;
	d__2 = *sa * d_imag(&cx[ix]);
	z__1.r = d__1, z__1.i = d__2;
	cx[i__2].r = z__1.r, cx[i__2].i = z__1.i;
	ix += *incx;
/* L10: */
    }
    return 0;

/*        CODE FOR INCREMENT EQUAL TO 1. */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	i__3 = i__;
	d__1 = *sa * cx[i__3].r;
	d__2 = *sa * d_imag(&cx[i__]);
	z__1.r = d__1, z__1.i = d__2;
	cx[i__2].r = z__1.r, cx[i__2].i = z__1.i;
/* L30: */
    }
    return 0;
} /* csscal_ */

/* -------------------------------------------------------------     ************ */
/*                                                                      CSUM */
/*                                                                  ************ */
/* Double Complex */ VOID csum_(doublecomplex * ret_val, integer *n, 
	doublecomplex *cx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, m, ix, mp1;
    static doublecomplex ctemp;


/*        TAKES THE SUM OF THE VALUES OF A VECTOR. */
/*        USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE. */
/*        JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --cx;

    /* Function Body */
     ret_val->r = 0.,  ret_val->i = 0.;
    ctemp.r = 0., ctemp.i = 0.;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ix;
	z__1.r = ctemp.r + cx[i__2].r, z__1.i = ctemp.i + cx[i__2].i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	ix += *incx;
/* L10: */
    }
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
    return ;

/*        CODE FOR INCREMENT EQUAL TO 1 */

/*        CLEAN-UP LOOP */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	z__1.r = ctemp.r + cx[i__2].r, z__1.i = ctemp.i + cx[i__2].i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 6) {
	i__2 = i__;
	z__6.r = ctemp.r + cx[i__2].r, z__6.i = ctemp.i + cx[i__2].i;
	i__3 = i__ + 1;
	z__5.r = z__6.r + cx[i__3].r, z__5.i = z__6.i + cx[i__3].i;
	i__4 = i__ + 2;
	z__4.r = z__5.r + cx[i__4].r, z__4.i = z__5.i + cx[i__4].i;
	i__5 = i__ + 3;
	z__3.r = z__4.r + cx[i__5].r, z__3.i = z__4.i + cx[i__5].i;
	i__6 = i__ + 4;
	z__2.r = z__3.r + cx[i__6].r, z__2.i = z__3.i + cx[i__6].i;
	i__7 = i__ + 5;
	z__1.r = z__2.r + cx[i__7].r, z__1.i = z__2.i + cx[i__7].i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
/* L50: */
    }
L60:
     ret_val->r = ctemp.r,  ret_val->i = ctemp.i;
    return ;
} /* csum_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     CSWAP */
/*                                                                  ************ */
/* Subroutine */ int cswap_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ctemp;


/*     INTERCHANGES TWO VECTORS. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*          CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS */
/*            NOT EQUAL TO 1. */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ix;
	ctemp.r = cx[i__2].r, ctemp.i = cx[i__2].i;
	i__2 = ix;
	i__3 = iy;
	cx[i__2].r = cy[i__3].r, cx[i__2].i = cy[i__3].i;
	i__2 = iy;
	cy[i__2].r = ctemp.r, cy[i__2].i = ctemp.i;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*          CODE FOR BOTH INCREMENTS EQUAL TO 1. */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	ctemp.r = cx[i__2].r, ctemp.i = cx[i__2].i;
	i__2 = i__;
	i__3 = i__;
	cx[i__2].r = cy[i__3].r, cx[i__2].i = cy[i__3].i;
	i__2 = i__;
	cy[i__2].r = ctemp.r, cy[i__2].i = ctemp.i;
/* L30: */
    }
    return 0;
} /* cswap_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     ICAMAX */
/*                                                                  ************ */
integer icamax_(integer *n, doublecomplex *cx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, ix;
    static doublereal smax;


/*     FINDS THE INDEX OF ELEMENT HAVING MAX. ABSOLUTE VALUE. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --cx;

    /* Function Body */
    ret_val = 0;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = ix;
    smax = (d__1 = cx[i__1].r, abs(d__1)) + (d__2 = d_imag(&cx[ix]), abs(d__2)
	    );
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = ix;
	if ((d__1 = cx[i__2].r, abs(d__1)) + (d__2 = d_imag(&cx[ix]), abs(
		d__2)) <= smax) {
	    goto L5;
	}
	ret_val = i__;
	i__2 = ix;
	smax = (d__1 = cx[i__2].r, abs(d__1)) + (d__2 = d_imag(&cx[ix]), abs(
		d__2));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENT EQUAL TO 1 */

L20:
    smax = (d__1 = cx[1].r, abs(d__1)) + (d__2 = d_imag(&cx[1]), abs(d__2));
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = i__;
	if ((d__1 = cx[i__2].r, abs(d__1)) + (d__2 = d_imag(&cx[i__]), abs(
		d__2)) <= smax) {
	    goto L30;
	}
	ret_val = i__;
	i__2 = i__;
	smax = (d__1 = cx[i__2].r, abs(d__1)) + (d__2 = d_imag(&cx[i__]), abs(
		d__2));
L30:
	;
    }
    return ret_val;
} /* icamax_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     ISAMIN */
/*                                                                  ************ */
integer isamin_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;
    doublereal d__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, ix;
    static doublereal dmin__;


/*        FINDS THE INDEX OF ELEMENT HAVING MIN. ABSOLUTE VALUE. */
/*        JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    dmin__ = (d__1 = dx[ix], abs(d__1));
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if ((d__1 = dx[ix], abs(d__1)) >= dmin__) {
	    goto L5;
	}
	ret_val = i__;
	dmin__ = (d__1 = dx[ix], abs(d__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENT EQUAL TO 1 */

L20:
    dmin__ = abs(dx[1]);
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if ((d__1 = dx[i__], abs(d__1)) >= dmin__) {
	    goto L30;
	}
	ret_val = i__;
	dmin__ = (d__1 = dx[i__], abs(d__1));
L30:
	;
    }
    return ret_val;
} /* isamin_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     ISMAX */
/*                                                                  ************ */
integer ismax_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, ix;
    static doublereal dmax__;


/*        FINDS THE INDEX OF ELEMENT WITH MAX. VALUE. */
/*        JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    dmax__ = dx[ix];
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (dx[ix] <= dmax__) {
	    goto L5;
	}
	ret_val = i__;
	dmax__ = dx[ix];
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENT EQUAL TO 1 */

L20:
    dmax__ = dx[1];
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (dx[i__] <= dmax__) {
	    goto L30;
	}
	ret_val = i__;
	dmax__ = dx[i__];
L30:
	;
    }
    return ret_val;
} /* ismax_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     ISMIN */
/*                                                                  ************ */
integer ismin_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, ix;
    static doublereal dmin__;


/*        FINDS THE INDEX OF ELEMENT WITH MIN. VALUE. */
/*        JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    dmin__ = dx[ix];
    ix += *incx;
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (dx[ix] >= dmin__) {
	    goto L5;
	}
	ret_val = i__;
	dmin__ = dx[ix];
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        CODE FOR INCREMENT EQUAL TO 1 */

L20:
    dmin__ = dx[1];
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (dx[i__] >= dmin__) {
	    goto L30;
	}
	ret_val = i__;
	dmin__ = dx[i__];
L30:
	;
    }
    return ret_val;
} /* ismin_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SCASUM */
/*                                                                  ************ */
doublereal scasum_(integer *n, doublecomplex *cx, integer *incx)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    double d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, ix;
    static doublereal stemp;


/*     TAKES THE SUM OF THE ABSOLUTE VALUES OF A COMPLEX VECTOR AND */
/*     RETURNS A SINGLE PRECISION RESULT. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --cx;

    /* Function Body */
    ret_val = 0.;
    stemp = 0.;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ix;
	stemp = stemp + (d__1 = cx[i__2].r, abs(d__1)) + (d__2 = d_imag(&cx[
		ix]), abs(d__2));
	ix += *incx;
/* L10: */
    }
    ret_val = stemp;
    return ret_val;

/*        CODE FOR INCREMENT EQUAL TO 1 */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	stemp = stemp + (d__1 = cx[i__2].r, abs(d__1)) + (d__2 = d_imag(&cx[
		i__]), abs(d__2));
/* L30: */
    }
    ret_val = stemp;
    return ret_val;
} /* scasum_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SCNRM2 */
/*                                                                  ************ */
doublereal scnrm2_(integer *n, doublecomplex *cx, integer *incx)
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal cutlo = 4.441e-16;
    static doublereal cuthi = 1.304e19;

    /* Format strings */
    static char fmt_30[] = "";
    static char fmt_50[] = "";
    static char fmt_70[] = "";
    static char fmt_90[] = "";
    static char fmt_110[] = "";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double d_imag(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, nn;
    static doublereal sum;
    static logical imag;
    static doublereal absx, xmax;
    static integer next;
    static logical scale;
    static doublereal hitest;

    /* Assigned format variables */
    static char *next_fmt;


    /* Parameter adjustments */
    --cx;

    /* Function Body */

/*     UNITARY NORM OF THE COMPLEX N-VECTOR STORED IN CX() WITH STORAGE */
/*     INCREMENT INCX . */
/*     IF    N .LE. 0 RETURN WITH RESULT = 0. */
/*     IF N .GE. 1 THEN INCX MUST BE .GE. 1 */

/*           C.L.LAWSON , 1978 JAN 08 */

/*     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE */
/*     HOPEFULLY APPLICABLE TO ALL MACHINES. */
/*         CUTLO = MAXIMUM OF  SQRT(U/EPS)  OVER ALL KNOWN MACHINES. */
/*         CUTHI = MINIMUM OF  SQRT(V)      OVER ALL KNOWN MACHINES. */
/*     WHERE */
/*         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1. */
/*         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT) */
/*         V   = LARGEST  NO.            (OVERFLOW  LIMIT) */

/*     BRIEF OUTLINE OF ALGORITHM.. */

/*     PHASE 1    SCANS ZERO COMPONENTS. */
/*     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO */
/*     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO */
/*     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M */
/*     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX. */

/*     VALUES FOR CUTLO AND CUTHI.. */
/*     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER */
/*     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS.. */
/*     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE */
/*                   UNIVAC AND DEC AT 2**(-103) */
/*                   THUS CUTLO = 2**(-51) = 4.44089E-16 */
/*     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC. */
/*                   THUS CUTHI = 2**(63.5) = 1.30438E19 */
/*     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC. */
/*                   THUS CUTLO = 2**(-33.5) = 8.23181D-11 */
/*     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19 */
/*     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 / */
/*     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 / */

    if (*n > 0) {
	goto L10;
    }
    ret_val = zero;
    goto L300;

L10:
    next = 0;
    next_fmt = fmt_30;
    sum = zero;
    nn = *n * *incx;
/*                    BEGIN MAIN LOOP */
    i__1 = nn;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	i__3 = i__;
	absx = (d__1 = cx[i__3].r, abs(d__1));
	imag = FALSE_;
	switch (next) {
	    case 0: goto L30;
	    case 1: goto L50;
	    case 2: goto L70;
	    case 3: goto L110;
	    case 4: goto L90;
	}
L30:
	if (absx > cutlo) {
	    goto L85;
	}
	next = 1;
	next_fmt = fmt_50;
	scale = FALSE_;

/*                        PHASE 1.  SUM IS ZERO */

L50:
	if (absx == zero) {
	    goto L200;
	}
	if (absx > cutlo) {
	    goto L85;
	}

/*                                PREPARE FOR PHASE 2. */
	next = 2;
	next_fmt = fmt_70;
	goto L105;

/*                                PREPARE FOR PHASE 4. */

L100:
	next = 3;
	next_fmt = fmt_110;
	sum = sum / absx / absx;
L105:
	scale = TRUE_;
	xmax = absx;
	goto L115;

/*                   PHASE 2.  SUM IS SMALL. */
/*                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW. */

L70:
	if (absx > cutlo) {
	    goto L75;
	}

/*                     COMMON CODE FOR PHASES 2 AND 4. */
/*                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW. */

L110:
	if (absx <= xmax) {
	    goto L115;
	}
/* Computing 2nd power */
	d__1 = xmax / absx;
	sum = one + sum * (d__1 * d__1);
	xmax = absx;
	goto L200;

L115:
/* Computing 2nd power */
	d__1 = absx / xmax;
	sum += d__1 * d__1;
	goto L200;


/*                  PREPARE FOR PHASE 3. */

L75:
	sum = sum * xmax * xmax;

L85:
	next = 4;
	next_fmt = fmt_90;
	scale = FALSE_;

/*     FOR REAL OR D.P. SET HITEST = CUTHI/N */
/*     FOR COMPLEX      SET HITEST = CUTHI/(2*N) */

	hitest = cuthi / (real) (*n);

/*                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING. */

L90:
	if (absx >= hitest) {
	    goto L100;
	}
/* Computing 2nd power */
	d__1 = absx;
	sum += d__1 * d__1;
L200:
/*                  CONTROL SELECTION OF REAL AND IMAGINARY PARTS. */

	if (imag) {
	    goto L210;
	}
	absx = (d__1 = d_imag(&cx[i__]), abs(d__1));
	imag = TRUE_;
	switch (next) {
	    case 0: goto L30;
	    case 1: goto L50;
	    case 2: goto L70;
	    case 3: goto L110;
	    case 4: goto L90;
	}

L210:
	;
    }

/*              END OF MAIN LOOP. */
/*              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING. */

    ret_val = sqrt(sum);
    if (scale) {
	ret_val *= xmax;
    }
L300:
    return ret_val;
} /* scnrm2_ */

/* -------------------------------------------------------------     ************ */
/*                                                                      SSUM */
/*                                                                  ************ */
doublereal ssum_(integer *n, doublereal *dx, integer *incx)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, m, ix, mp1;
    static doublereal dtemp;


/*        TAKES THE SUM OF THE VALUES OF A VECTOR. */
/*        USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE. */
/*        JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[ix];
	ix += *incx;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        CODE FOR INCREMENT EQUAL TO 1 */

/*        CLEAN-UP LOOP */

L20:
    m = *n % 6;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[i__];
/* L30: */
    }
    if (*n < 6) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 6) {
	dtemp = dtemp + dx[i__] + dx[i__ + 1] + dx[i__ + 2] + dx[i__ + 3] + 
		dx[i__ + 4] + dx[i__ + 5];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ssum_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SROTMG */
/*                                                                  ************ */
/* Subroutine */ int srotmg_(doublereal *d1, doublereal *d2, doublereal *b1, 
	doublereal *b2, doublereal *param)
{
    /* System generated locals */
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal h__[4]	/* was [2][2] */;
    static integer i__;
    static doublereal m, u, temp;
    extern /* Subroutine */ int srotmgw_(doublereal *, integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *);


/*     This subroutine computes a modified Givens plane rotation matrix. */
/*     Martin J. McBride, 7/10/85. */
/*     General Electric CRD, Information System Operation. */


    /* Parameter adjustments */
    --param;

    /* Function Body */
    m = 4096.;
/*  Case 4:  If D2 or B2 are equal to 0, then PARAM(1) = 2 is the only */
/*           change. */
    if (*d2 * *b2 == 0.) {
	param[1] = -2.;
	return 0;
    }
/*  If D1 is less than 0, then PARAM(1) = -1 and PARAM(2) to PARAM(5) */
/*           are set to 0. */
    if (*d1 < 0.) {
	*d1 = 0.;
	*d2 = 0.;
	*b1 = 0.;
	param[1] = -1.;
	for (i__ = 2; i__ <= 5; ++i__) {
	    param[i__] = 0.;
/* L5: */
	}
	return 0;
    }
    if ((d__1 = *d1 * *b1 * *b1, abs(d__1)) > (d__2 = *d2 * *b2 * *b2, abs(
	    d__2))) {
/*  Case 1:  D1 and D2 are updated by a factor of 1/U, where U is the */
/*           determinant of matrix H. */
	h__[0] = 1.;
	h__[2] = *d2 * *b2 / (*d1 * *b1);
	h__[1] = -(*b2) / *b1;
	h__[3] = 1.;
	param[1] = 0.;
	param[3] = h__[1];
	param[4] = h__[2];
	u = *d2 * *b2 * *b2 / (*d1 * *b1 * *b1) + 1.;
	*d1 /= u;
	*d2 /= u;
	*b1 *= u;
	if (*d1 == 0. || *d2 == 0.) {
	    return 0;
	}
	if (abs(*d1) < 1.f / (m * m) || abs(*d1) > m * m) {
	    srotmgw_(d1, &c__1, b1, h__, &param[1], &m);
	}
	if (abs(*d2) < 1.f / (m * m) || abs(*d2) > m * m) {
	    srotmgw_(d2, &c__2, b1, h__, &param[1], &m);
	}
    } else {
/*    If D2 is less than 0, then PARAM(1) = -1 and PARAM(2) to PARAM(5) */
/*             are set to 0. */
	if (*d2 < 0.) {
	    *d1 = 0.;
	    *d2 = 0.;
	    *b1 = 0.;
	    param[1] = -1.;
	    for (i__ = 2; i__ <= 5; ++i__) {
		param[i__] = 0.;
/* L15: */
	    }
	    return 0;
	}
/*  Case 2:  D1 and D2 are updated by a factor of 1/U, where U is the */
/*           determinant of matrix H, and then are interchanged. */
	h__[0] = *d1 * *b1 / (*d2 * *b2);
	h__[2] = 1.;
	h__[1] = -1.;
	h__[3] = *b1 / *b2;
	param[1] = 1.;
	param[2] = h__[0];
	param[5] = h__[3];
	u = *d1 * *b1 * *b1 / (*d2 * *b2 * *b2) + 1.;
	temp = *d1 / u;
	*d1 = *d2 / u;
	*d2 = temp;
	*b1 = *b2 * u;
	if (*d1 == 0. || *d2 == 0.) {
	    return 0;
	}
	if (abs(*d1) < 1.f / (m * m) || abs(*d1) > m * m) {
	    srotmgw_(d1, &c__1, b1, h__, &param[1], &m);
	}
	if (abs(*d2) < 1.f / (m * m) || abs(*d2) > m * m) {
	    srotmgw_(d2, &c__2, b1, h__, &param[1], &m);
	}
    }
    return 0;
} /* srotmg_ */

/* Subroutine */ int srotmgw_(doublereal *di, integer *i__, doublereal *b1, 
	doublereal *h__, doublereal *param, doublereal *m)
{
/*  Case 3:  D1 or D2 is updated until it is within a given window, */
/*           m**-2 <= abs(D(i)) <= m**2 where m=4096. */
    /* Parameter adjustments */
    --param;
    h__ -= 3;

    /* Function Body */
L10:
    if (abs(*di) < 1.f / (*m * *m)) {
	*di = *di * *m * *m;
	if (*i__ == 1) {
	    *b1 /= *m;
	}
	h__[*i__ + 2] /= *m;
	h__[*i__ + 4] /= *m;
    } else if (abs(*di) > *m * *m) {
	*di = *di * 1.f / (*m * *m);
	if (*i__ == 1) {
	    *b1 *= *m;
	}
	h__[*i__ + 2] *= *m;
	h__[*i__ + 4] *= *m;
    }
    if (abs(*di) < 1.f / (*m * *m) || abs(*di) > *m * *m) {
	goto L10;
    }
    param[1] = -1.;
    param[2] = h__[3];
    param[3] = h__[4];
    param[4] = h__[5];
    param[5] = h__[6];
    return 0;
} /* srotmgw_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     SROTM */
/*                                                                  ************ */
/* Subroutine */ int srotm_(integer *n, doublereal *sx, integer *incx, 
	doublereal *sy, integer *incy, doublereal *param)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsle(cilist *), e_wsle(void), do_lio(integer *, integer *, char 
	    *, ftnlen);

    /* Local variables */
    static doublereal h__[4]	/* was [2][2] */;
    static integer i__;
    static real x, y;
    static integer ix, iy;

    /* Fortran I/O blocks */
    static cilist io___138 = { 0, 6, 0, 0, 0 };
    static cilist io___139 = { 0, 6, 0, 0, 0 };
    static cilist io___140 = { 0, 6, 0, 0, 0 };



/*     This subroutine applies the modified Givens rotation matrix. */
/*     Martin J. McBride.  7/11/85. */
/*     General Electric CRD, Information System Operation. */

    /* Parameter adjustments */
    --param;
    --sy;
    --sx;

    /* Function Body */
    if (*n < 0) {
	s_stop("", (ftnlen)0);
    }
    if (*n == 0) {
	return 0;
    }
    if (param[1] == -2.) {
	return 0;
    }
/*  Conditions for setting up matrix H for multiplication. */
    if (param[1] == 1.) {
	h__[0] = param[2];
	h__[1] = -1.;
	h__[2] = 1.;
	h__[3] = param[5];
    } else if (param[1] == 0.) {
	h__[0] = 1.;
	h__[1] = param[3];
	h__[2] = param[4];
	h__[3] = 1.;
    } else if (param[1] == -1.) {
	h__[0] = param[2];
	h__[1] = param[3];
	h__[2] = param[4];
	h__[3] = param[5];
    } else {
	s_wsle(&io___138);
	e_wsle();
	s_wsle(&io___139);
	do_lio(&c__9, &c__1, "     SROTM called with incorrect parameter key",
		 (ftnlen)46);
	e_wsle();
	s_wsle(&io___140);
	e_wsle();
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        Code for unequal increments of vectors X and Y or */
/*          equal increments not equal to 1. */

    ix = 1;
    iy = 1;
    if (*incx <= 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy <= 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x = sx[ix];
	y = sy[iy];
	sx[ix] = x * h__[0] + y * h__[2];
	sy[iy] = x * h__[1] + y * h__[3];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        Code for equal increments of vectors X and Y. */

L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x = sx[i__];
	y = sy[i__];
	sx[i__] = x * h__[0] + y * h__[2];
	sy[i__] = x * h__[1] + y * h__[3];
/* L30: */
    }
    return 0;
} /* srotm_ */

