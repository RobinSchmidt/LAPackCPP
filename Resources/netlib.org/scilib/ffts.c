/* ffts.f -- translated by f2c (version 20100827).
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

static integer c__9 = 9;
static integer c__1 = 1;

/* -------------------------------------------------------------     ************ */
/*                                                                     CRFFT2 */
/*                                                                  ************ */
/* Subroutine */ int crfft2_(integer *init, integer *ix, integer *n, complex *
	cx, complex *cwork, real *cy)
{
    static integer l2, n2, nn, ls, ns;
    extern /* Subroutine */ int abort_(char *, ftnlen), crble1_(integer *, 
	    complex *), crock1_(integer *, complex *, complex *), crock2_(), 
	    crock3_(), crform_(integer *, integer *, integer *, complex *, 
	    complex *, complex *);
    extern integer popcnt_(integer *);


/*     THE STOCKHAM AUTO-SORT FFT */


    /* Parameter adjustments */
    --cy;
    --cwork;
    --cx;

    /* Function Body */
    l2 = popcnt_(n);
    if (l2 != 1) {
	abort_("CRFFT2", (ftnlen)6);
    }
    n2 = *n + 3;
    nn = *n / 2;
    if (*init == 0) {
	goto L10;
    }
    if (*n < 4) {
	abort_("CRFFT2", (ftnlen)6);
    }
    crble1_(&nn, &cwork[n2]);
    return 0;
L10:
    ns = *n / 4;
    crform_(ix, &ns, &nn, &cx[1], &cwork[nn + 1], &cwork[n2]);
    crock1_(&ns, &cwork[1], &cwork[nn + 1]);
    if (*ix <= 0) {
	goto L50;
    }
    ls = 2;
    ns /= 2;
L20:
    if (ns == 1) {
	goto L30;
    }
    crock2_(&ls, &ns, &cwork[nn + 1], &cwork[1], &cwork[n2]);
    ls += ls;
    ns /= 2;
    if (ns == 1) {
	goto L130;
    }
    crock2_(&ls, &ns, &cwork[1], &cwork[nn + 1], &cwork[n2]);
    ls += ls;
    ns /= 2;
    goto L20;
L30:
    crock2_(&ls, &ns, &cy[1], &cwork[1], &cwork[n2]);
    return 0;
L130:
    crock2_(&ls, &ns, &cy[1], &cwork[nn + 1], &cwork[n2]);
    return 0;
L50:
    ls = 2;
    ns /= 2;
L60:
    if (ns == 1) {
	goto L70;
    }
    crock3_(&ls, &ns, &cwork[nn + 1], &cwork[1], &cwork[n2]);
    ls += ls;
    ns /= 2;
    if (ns == 1) {
	goto L120;
    }
    crock3_(&ls, &ns, &cwork[1], &cwork[nn + 1], &cwork[n2]);
    ls += ls;
    ns /= 2;
    goto L60;
L70:
    crock3_(&ls, &ns, &cy[1], &cwork[1], &cwork[n2]);
    return 0;
L120:
    crock3_(&ls, &ns, &cy[1], &cwork[nn + 1], &cwork[n2]);
    return 0;
} /* crfft2_ */

/* -----------------------------------------------  ************ */
/*                                                    CRFORM */
/*                                                 ************ */
/* Subroutine */ int crform_(integer *ix, integer *ns, integer *ndiv2, 
	complex *cx, complex *c__, real *ch2)
{
    /* System generated locals */
    integer c_dim1, c_offset, ch2_dim1, ch2_offset, i__1, i__2, i__3, i__4, 
	    i__5, i__6;
    real r__1;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer i__, k;
    static complex wyk, wyk1;



    /* Parameter adjustments */
    c_dim1 = *ns;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;
    ch2_dim1 = *ndiv2;
    ch2_offset = 1 + ch2_dim1;
    ch2 -= ch2_offset;
    --cx;

    /* Function Body */
    if (*ix > 0) {
	goto L50;
    }
    k = *ns + 1;
    i__1 = *ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r_cnjg(&q__1, &cx[*ndiv2 - i__ + 2]);
	wyk.r = q__1.r, wyk.i = q__1.i;
	i__2 = i__ + c_dim1;
	i__3 = i__;
	q__2.r = cx[i__3].r + wyk.r, q__2.i = cx[i__3].i + wyk.i;
	i__4 = i__;
	q__4.r = cx[i__4].r - wyk.r, q__4.i = cx[i__4].i - wyk.i;
	i__5 = i__ + (ch2_dim1 << 1);
	i__6 = i__ + ch2_dim1;
	q__5.r = ch2[i__5], q__5.i = ch2[i__6];
	q__3.r = q__4.r * q__5.r - q__4.i * q__5.i, q__3.i = q__4.r * q__5.i 
		+ q__4.i * q__5.r;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
	r_cnjg(&q__1, &cx[*ndiv2 - k + 2]);
	wyk1.r = q__1.r, wyk1.i = q__1.i;
	i__2 = i__ + (c_dim1 << 1);
	i__3 = k;
	q__2.r = cx[i__3].r + wyk1.r, q__2.i = cx[i__3].i + wyk1.i;
	i__4 = k;
	q__4.r = cx[i__4].r - wyk1.r, q__4.i = cx[i__4].i - wyk1.i;
	i__5 = k + (ch2_dim1 << 1);
	i__6 = k + ch2_dim1;
	q__5.r = ch2[i__5], q__5.i = ch2[i__6];
	q__3.r = q__4.r * q__5.r - q__4.i * q__5.i, q__3.i = q__4.r * q__5.i 
		+ q__4.i * q__5.r;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
	++k;
/* L10: */
    }
    return 0;
L50:
    k = *ns + 1;
    i__1 = *ns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r_cnjg(&q__1, &cx[*ndiv2 - i__ + 2]);
	wyk.r = q__1.r, wyk.i = q__1.i;
	i__2 = i__ + c_dim1;
	i__3 = i__;
	q__2.r = cx[i__3].r + wyk.r, q__2.i = cx[i__3].i + wyk.i;
	i__4 = i__;
	q__4.r = cx[i__4].r - wyk.r, q__4.i = cx[i__4].i - wyk.i;
	r__1 = -ch2[i__ + (ch2_dim1 << 1)];
	i__5 = i__ + ch2_dim1;
	q__5.r = r__1, q__5.i = ch2[i__5];
	q__3.r = q__4.r * q__5.r - q__4.i * q__5.i, q__3.i = q__4.r * q__5.i 
		+ q__4.i * q__5.r;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
	r_cnjg(&q__1, &cx[*ndiv2 - k + 2]);
	wyk1.r = q__1.r, wyk1.i = q__1.i;
	i__2 = i__ + (c_dim1 << 1);
	i__3 = k;
	q__2.r = cx[i__3].r + wyk1.r, q__2.i = cx[i__3].i + wyk1.i;
	i__4 = k;
	q__4.r = cx[i__4].r - wyk1.r, q__4.i = cx[i__4].i - wyk1.i;
	r__1 = -ch2[k + (ch2_dim1 << 1)];
	i__5 = k + ch2_dim1;
	q__5.r = r__1, q__5.i = ch2[i__5];
	q__3.r = q__4.r * q__5.r - q__4.i * q__5.i, q__3.i = q__4.r * q__5.i 
		+ q__4.i * q__5.r;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
	++k;
/* L20: */
    }
    return 0;
} /* crform_ */

/* -----------------------------------------------  ************ */
/*                                                    CROCK1 */
/*                                                 ************ */
/* Subroutine */ int crock1_(integer *ns, complex *c__, complex *ch)
{
    /* System generated locals */
    integer c_dim1, c_offset, ch_dim1, ch_offset, i__1, i__2, i__3, i__4;
    complex q__1;

    /* Local variables */
    static integer j;



    /* Parameter adjustments */
    ch_dim1 = *ns;
    ch_offset = 1 + ch_dim1;
    ch -= ch_offset;
    c_dim1 = *ns;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    i__1 = *ns;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + c_dim1;
	i__3 = j + ch_dim1;
	i__4 = j + (ch_dim1 << 1);
	q__1.r = ch[i__3].r + ch[i__4].r, q__1.i = ch[i__3].i + ch[i__4].i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
	i__2 = j + (c_dim1 << 1);
	i__3 = j + ch_dim1;
	i__4 = j + (ch_dim1 << 1);
	q__1.r = ch[i__3].r - ch[i__4].r, q__1.i = ch[i__3].i - ch[i__4].i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
/* L300: */
    }
    return 0;
} /* crock1_ */

/* -----------------------------------------------  ************ */
/*                                                    CROCK2 */
/*                                                 ************ */
/* Subroutine */ int crock2_(integer *ls, integer *ns, complex *c__, complex *
	ch, real *ch2)
{
    /* System generated locals */
    integer c_dim1, c_dim2, c_offset, ch_dim1, ch_offset, ch2_dim2, ch2_dim3, 
	    ch2_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2;

    /* Local variables */
    static integer i__, j;
    static complex wyk;



    /* Parameter adjustments */
    ch2_dim2 = *ns;
    ch2_dim3 = *ls;
    ch2_offset = 1 + 2 * (1 + ch2_dim2 * (1 + ch2_dim3));
    ch2 -= ch2_offset;
    ch_dim1 = *ns;
    ch_offset = 1 + ch_dim1 * 3;
    ch -= ch_offset;
    c_dim1 = *ns;
    c_dim2 = *ls;
    c_offset = 1 + c_dim1 * (1 + c_dim2);
    c__ -= c_offset;

    /* Function Body */
    if (*ls > *ns) {
	goto L20;
    }
    i__1 = *ls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ns;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = ((i__ + ch2_dim3) * ch2_dim2 + 1 << 1) + 1;
	    i__4 = ((i__ + (ch2_dim3 << 1)) * ch2_dim2 + 1 << 1) + 1;
	    q__2.r = ch2[i__3], q__2.i = ch2[i__4];
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__1.r = q__2.r * ch[i__5].r - q__2.i * ch[i__5].i, q__1.i = 
		    q__2.r * ch[i__5].i + q__2.i * ch[i__5].r;
	    wyk.r = q__1.r, wyk.i = q__1.i;
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r + wyk.r, q__1.i = ch[i__4].i + wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r - wyk.r, q__1.i = ch[i__4].i - wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L200: */
	}
    }
    return 0;
L20:
    i__2 = *ns;
    for (j = 1; j <= i__2; ++j) {
	i__1 = *ls;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = ((i__ + ch2_dim3) * ch2_dim2 + 1 << 1) + 1;
	    i__4 = ((i__ + (ch2_dim3 << 1)) * ch2_dim2 + 1 << 1) + 1;
	    q__2.r = ch2[i__3], q__2.i = ch2[i__4];
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__1.r = q__2.r * ch[i__5].r - q__2.i * ch[i__5].i, q__1.i = 
		    q__2.r * ch[i__5].i + q__2.i * ch[i__5].r;
	    wyk.r = q__1.r, wyk.i = q__1.i;
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r + wyk.r, q__1.i = ch[i__4].i + wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r - wyk.r, q__1.i = ch[i__4].i - wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L400: */
	}
    }
    return 0;
} /* crock2_ */

/* -----------------------------------------------  ************ */
/*                                                    CROCK3 */
/*                                                 ************ */
/* Subroutine */ int crock3_(integer *ls, integer *ns, complex *c__, complex *
	ch, real *ch2)
{
    /* System generated locals */
    integer c_dim1, c_dim2, c_offset, ch_dim1, ch_offset, ch2_dim2, ch2_dim3, 
	    ch2_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer i__, j;
    static complex wyk;



    /* Parameter adjustments */
    ch2_dim2 = *ns;
    ch2_dim3 = *ls;
    ch2_offset = 1 + 2 * (1 + ch2_dim2 * (1 + ch2_dim3));
    ch2 -= ch2_offset;
    ch_dim1 = *ns;
    ch_offset = 1 + ch_dim1 * 3;
    ch -= ch_offset;
    c_dim1 = *ns;
    c_dim2 = *ls;
    c_offset = 1 + c_dim1 * (1 + c_dim2);
    c__ -= c_offset;

    /* Function Body */
    if (*ls > *ns) {
	goto L30;
    }
    i__1 = *ls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ns;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = ((i__ + ch2_dim3) * ch2_dim2 + 1 << 1) + 1;
	    i__4 = ((i__ + (ch2_dim3 << 1)) * ch2_dim2 + 1 << 1) + 1;
	    q__3.r = ch2[i__3], q__3.i = ch2[i__4];
	    r_cnjg(&q__2, &q__3);
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__1.r = q__2.r * ch[i__5].r - q__2.i * ch[i__5].i, q__1.i = 
		    q__2.r * ch[i__5].i + q__2.i * ch[i__5].r;
	    wyk.r = q__1.r, wyk.i = q__1.i;
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r + wyk.r, q__1.i = ch[i__4].i + wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r - wyk.r, q__1.i = ch[i__4].i - wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L600: */
	}
    }
    return 0;
L30:
    i__2 = *ns;
    for (j = 1; j <= i__2; ++j) {
	i__1 = *ls;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = ((i__ + ch2_dim3) * ch2_dim2 + 1 << 1) + 1;
	    i__4 = ((i__ + (ch2_dim3 << 1)) * ch2_dim2 + 1 << 1) + 1;
	    q__3.r = ch2[i__3], q__3.i = ch2[i__4];
	    r_cnjg(&q__2, &q__3);
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__1.r = q__2.r * ch[i__5].r - q__2.i * ch[i__5].i, q__1.i = 
		    q__2.r * ch[i__5].i + q__2.i * ch[i__5].r;
	    wyk.r = q__1.r, wyk.i = q__1.i;
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r + wyk.r, q__1.i = ch[i__4].i + wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r - wyk.r, q__1.i = ch[i__4].i - wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L800: */
	}
    }
    return 0;
} /* crock3_ */

/* -----------------------------------------------  ************ */
/*                                                    CRBLE1 */
/*                                                 ************ */
/* Subroutine */ int crble1_(integer *nn, real *work)
{
    /* Initialized data */

    static real twopi = 6.28318530717958647692f;

    /* System generated locals */
    integer work_dim1, work_offset, i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, n;
    static real p2;


    /* Parameter adjustments */
    work_dim1 = *nn;
    work_offset = 1 + work_dim1;
    work -= work_offset;

    /* Function Body */

    n = 2 * *nn;
    p2 = twopi / n;
    i__1 = *nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[i__ + work_dim1] = cos(p2 * (i__ - 1));
	work[i__ + (work_dim1 << 1)] = sin(p2 * (i__ - 1));
/* L10: */
    }
    return 0;
} /* crble1_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     RCFFT2 */
/*                                                                  ************ */
/* Subroutine */ int rcfft2_(integer *init, integer *ix, integer *n, complex *
	cx, complex *cwork, complex *cy)
{
    static integer l2, n2, nn, ls, ns;
    extern /* Subroutine */ int abort_(char *, ftnlen), rable1_(integer *, 
	    complex *), rtock1_(integer *, complex *, complex *), rtock2_(
	    integer *, integer *, complex *, complex *, complex *), rtock3_(
	    integer *, integer *, complex *, complex *, complex *), rconv1_(
	    integer *, complex *, complex *, complex *), rconv2_(integer *, 
	    complex *, complex *, complex *);
    extern integer popcnt_(integer *);


/*     THE STOCKHAM AUTO-SORT FFT */


    /* Parameter adjustments */
    --cy;
    --cwork;
    --cx;

    /* Function Body */
    l2 = popcnt_(n);
    if (l2 != 1) {
	abort_("RCFFT2", (ftnlen)6);
    }
    n2 = *n + 3;
    nn = *n / 2;
    if (*init == 0) {
	goto L10;
    }
    if (*n < 4) {
	abort_("RCFFT2", (ftnlen)6);
    }
    rable1_(&nn, &cwork[n2]);
    return 0;
L10:
    ns = *n / 4;
    rtock1_(&ns, &cwork[1], &cx[1]);
    if (*ix <= 0) {
	goto L50;
    }
    ls = 2;
    ns /= 2;
L20:
    rtock2_(&ls, &ns, &cwork[nn + 1], &cwork[1], &cwork[n2]);
    if (ns == 1) {
	goto L30;
    }
    ls += ls;
    ns /= 2;
    rtock2_(&ls, &ns, &cwork[1], &cwork[nn + 1], &cwork[n2]);
    if (ns == 1) {
	goto L130;
    }
    ls += ls;
    ns /= 2;
    goto L20;
L30:
    rconv1_(n, &cy[1], &cwork[nn + 1], &cwork[n2]);
    return 0;
L130:
    rconv1_(n, &cy[1], &cwork[1], &cwork[n2]);
    return 0;
L50:
    ls = 2;
    ns /= 2;
L60:
    rtock3_(&ls, &ns, &cwork[nn + 1], &cwork[1], &cwork[n2]);
    if (ns == 1) {
	goto L70;
    }
    ls += ls;
    ns /= 2;
    rtock3_(&ls, &ns, &cwork[1], &cwork[nn + 1], &cwork[n2]);
    if (ns == 1) {
	goto L120;
    }
    ls += ls;
    ns /= 2;
    goto L60;
L70:
    rconv2_(n, &cy[1], &cwork[nn + 1], &cwork[n2]);
    return 0;
L120:
    rconv2_(n, &cy[1], &cwork[1], &cwork[n2]);
    return 0;
} /* rcfft2_ */

/* -----------------------------------------------  ************ */
/*                                                    RTOCK2 */
/*                                                 ************ */
/* Subroutine */ int rtock2_(integer *ls, integer *ns, complex *c__, complex *
	ch, real *ch2)
{
    /* System generated locals */
    integer c_dim1, c_dim2, c_offset, ch_dim1, ch_offset, ch2_dim2, ch2_dim3, 
	    ch2_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2;

    /* Local variables */
    static integer i__, j;
    static complex wyk;



    /* Parameter adjustments */
    ch2_dim2 = *ns;
    ch2_dim3 = *ls;
    ch2_offset = 1 + 2 * (1 + ch2_dim2 * (1 + ch2_dim3));
    ch2 -= ch2_offset;
    ch_dim1 = *ns;
    ch_offset = 1 + ch_dim1 * 3;
    ch -= ch_offset;
    c_dim1 = *ns;
    c_dim2 = *ls;
    c_offset = 1 + c_dim1 * (1 + c_dim2);
    c__ -= c_offset;

    /* Function Body */
    if (*ls > *ns) {
	goto L20;
    }
    i__1 = *ls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ns;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = ((i__ + ch2_dim3) * ch2_dim2 + 1 << 1) + 1;
	    i__4 = ((i__ + (ch2_dim3 << 1)) * ch2_dim2 + 1 << 1) + 1;
	    q__2.r = ch2[i__3], q__2.i = ch2[i__4];
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__1.r = q__2.r * ch[i__5].r - q__2.i * ch[i__5].i, q__1.i = 
		    q__2.r * ch[i__5].i + q__2.i * ch[i__5].r;
	    wyk.r = q__1.r, wyk.i = q__1.i;
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r + wyk.r, q__1.i = ch[i__4].i + wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r - wyk.r, q__1.i = ch[i__4].i - wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L200: */
	}
    }
    return 0;
L20:
    i__2 = *ns;
    for (j = 1; j <= i__2; ++j) {
	i__1 = *ls;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = ((i__ + ch2_dim3) * ch2_dim2 + 1 << 1) + 1;
	    i__4 = ((i__ + (ch2_dim3 << 1)) * ch2_dim2 + 1 << 1) + 1;
	    q__2.r = ch2[i__3], q__2.i = ch2[i__4];
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__1.r = q__2.r * ch[i__5].r - q__2.i * ch[i__5].i, q__1.i = 
		    q__2.r * ch[i__5].i + q__2.i * ch[i__5].r;
	    wyk.r = q__1.r, wyk.i = q__1.i;
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r + wyk.r, q__1.i = ch[i__4].i + wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r - wyk.r, q__1.i = ch[i__4].i - wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L400: */
	}
    }
    return 0;
} /* rtock2_ */

/* -----------------------------------------------  ************ */
/*                                                    RTOCK1 */
/*                                                 ************ */
/* Subroutine */ int rtock1_(integer *ns, complex *c__, complex *ch)
{
    /* System generated locals */
    integer c_dim1, c_offset, ch_dim1, ch_offset, i__1, i__2, i__3, i__4;
    complex q__1;

    /* Local variables */
    static integer j;



    /* Parameter adjustments */
    ch_dim1 = *ns;
    ch_offset = 1 + ch_dim1 * 3;
    ch -= ch_offset;
    c_dim1 = *ns;
    c_offset = 1 + (c_dim1 << 1);
    c__ -= c_offset;

    /* Function Body */
    i__1 = *ns;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + (c_dim1 << 1);
	i__3 = j + ch_dim1 * 3;
	i__4 = j + (ch_dim1 << 2);
	q__1.r = ch[i__3].r + ch[i__4].r, q__1.i = ch[i__3].i + ch[i__4].i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
	i__2 = j + c_dim1 * 3;
	i__3 = j + ch_dim1 * 3;
	i__4 = j + (ch_dim1 << 2);
	q__1.r = ch[i__3].r - ch[i__4].r, q__1.i = ch[i__3].i - ch[i__4].i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
/* L300: */
    }
    return 0;
} /* rtock1_ */

/* -----------------------------------------------  ************ */
/*                                                    RTOCK3 */
/*                                                 ************ */
/* Subroutine */ int rtock3_(integer *ls, integer *ns, complex *c__, complex *
	ch, real *ch2)
{
    /* System generated locals */
    integer c_dim1, c_dim2, c_offset, ch_dim1, ch_offset, ch2_dim2, ch2_dim3, 
	    ch2_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer i__, j;
    static complex wyk;



    /* Parameter adjustments */
    ch2_dim2 = *ns;
    ch2_dim3 = *ls;
    ch2_offset = 1 + 2 * (1 + ch2_dim2 * (1 + ch2_dim3));
    ch2 -= ch2_offset;
    ch_dim1 = *ns;
    ch_offset = 1 + ch_dim1 * 3;
    ch -= ch_offset;
    c_dim1 = *ns;
    c_dim2 = *ls;
    c_offset = 1 + c_dim1 * (1 + c_dim2);
    c__ -= c_offset;

    /* Function Body */
    if (*ls > *ns) {
	goto L30;
    }
    i__1 = *ls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ns;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = ((i__ + ch2_dim3) * ch2_dim2 + 1 << 1) + 1;
	    i__4 = ((i__ + (ch2_dim3 << 1)) * ch2_dim2 + 1 << 1) + 1;
	    q__3.r = ch2[i__3], q__3.i = ch2[i__4];
	    r_cnjg(&q__2, &q__3);
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__1.r = q__2.r * ch[i__5].r - q__2.i * ch[i__5].i, q__1.i = 
		    q__2.r * ch[i__5].i + q__2.i * ch[i__5].r;
	    wyk.r = q__1.r, wyk.i = q__1.i;
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r + wyk.r, q__1.i = ch[i__4].i + wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r - wyk.r, q__1.i = ch[i__4].i - wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L600: */
	}
    }
    return 0;
L30:
    i__2 = *ns;
    for (j = 1; j <= i__2; ++j) {
	i__1 = *ls;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = ((i__ + ch2_dim3) * ch2_dim2 + 1 << 1) + 1;
	    i__4 = ((i__ + (ch2_dim3 << 1)) * ch2_dim2 + 1 << 1) + 1;
	    q__3.r = ch2[i__3], q__3.i = ch2[i__4];
	    r_cnjg(&q__2, &q__3);
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__1.r = q__2.r * ch[i__5].r - q__2.i * ch[i__5].i, q__1.i = 
		    q__2.r * ch[i__5].i + q__2.i * ch[i__5].r;
	    wyk.r = q__1.r, wyk.i = q__1.i;
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r + wyk.r, q__1.i = ch[i__4].i + wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    q__1.r = ch[i__4].r - wyk.r, q__1.i = ch[i__4].i - wyk.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L800: */
	}
    }
    return 0;
} /* rtock3_ */

/* -----------------------------------------------  ************ */
/*                                                    CRFORM */
/*                                                 ************ */
/* Subroutine */ int rable1_(integer *nn, real *work)
{
    /* Initialized data */

    static real twopi = 6.28318530717958647692f;

    /* System generated locals */
    integer work_dim1, work_offset, i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, n;
    static real p2;


    /* Parameter adjustments */
    work_dim1 = *nn;
    work_offset = 1 + work_dim1;
    work -= work_offset;

    /* Function Body */

    n = 2 * *nn;
    p2 = twopi / n;
    i__1 = *nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[i__ + work_dim1] = cos(p2 * (i__ - 1));
	work[i__ + (work_dim1 << 1)] = sin(p2 * (i__ - 1));
/* L10: */
    }
    return 0;
} /* rable1_ */

/* -----------------------------------------------  ************ */
/*                                                    RCONV1 */
/*                                                 ************ */
/* Subroutine */ int rconv1_(integer *n, complex *cy, real *c__, real *ch)
{
    /* System generated locals */
    integer ch_dim1, ch_offset, i__1, i__2;
    complex q__1;

    /* Local variables */
    static integer i__, k;
    static real p[2]	/* was [2][1] */, x, y, z__;
    static integer n2;
    static real z1;



    /* Parameter adjustments */
    ch_dim1 = *n / 2;
    ch_offset = 1 + ch_dim1;
    ch -= ch_offset;
    c__ -= 3;
    --cy;

    /* Function Body */
    n2 = *n / 2;
    p[0] = (c__[3] + c__[4]) * 2;
    p[1] = (c__[3] - c__[4]) * 2;
    q__1.r = p[0], q__1.i = 0.f;
    cy[1].r = q__1.r, cy[1].i = q__1.i;
    i__1 = n2 + 1;
    q__1.r = p[1], q__1.i = 0.f;
    cy[i__1].r = q__1.r, cy[i__1].i = q__1.i;
    k = n2;
    i__1 = n2;
    for (i__ = 2; i__ <= i__1; ++i__) {
	x = c__[(i__ << 1) + 1] + c__[(k << 1) + 1];
	y = c__[(i__ << 1) + 2] + c__[(k << 1) + 2];
	z__ = c__[(i__ << 1) + 1] - c__[(k << 1) + 1];
	z1 = c__[(i__ << 1) + 2] - c__[(k << 1) + 2];
	p[0] = x + ch[i__ + ch_dim1] * y + ch[i__ + (ch_dim1 << 1)] * z__;
	p[1] = z1 + ch[i__ + (ch_dim1 << 1)] * y - ch[i__ + ch_dim1] * z__;
	i__2 = i__;
	q__1.r = p[0], q__1.i = p[1];
	cy[i__2].r = q__1.r, cy[i__2].i = q__1.i;
	--k;
/* L10: */
    }
    return 0;
} /* rconv1_ */

/* -----------------------------------------------  ************ */
/*                                                    RCONV2 */
/*                                                 ************ */
/* Subroutine */ int rconv2_(integer *n, complex *cy, real *c__, real *ch)
{
    /* System generated locals */
    integer ch_dim1, ch_offset, i__1, i__2;
    complex q__1;

    /* Local variables */
    static integer i__, k;
    static real p[2]	/* was [2][1] */, x, y, z__;
    static integer n2;
    static real z1;



    /* Parameter adjustments */
    ch_dim1 = *n / 2;
    ch_offset = 1 + ch_dim1;
    ch -= ch_offset;
    c__ -= 3;
    --cy;

    /* Function Body */
    n2 = *n / 2;
    p[0] = (c__[3] + c__[4]) * 2;
    p[1] = (c__[3] - c__[4]) * 2;
    q__1.r = p[0], q__1.i = 0.f;
    cy[1].r = q__1.r, cy[1].i = q__1.i;
    i__1 = n2 + 1;
    q__1.r = p[1], q__1.i = 0.f;
    cy[i__1].r = q__1.r, cy[i__1].i = q__1.i;
    k = n2;
    i__1 = n2;
    for (i__ = 2; i__ <= i__1; ++i__) {
	x = c__[(i__ << 1) + 1] + c__[(k << 1) + 1];
	y = c__[(i__ << 1) + 2] + c__[(k << 1) + 2];
	z__ = c__[(i__ << 1) + 1] - c__[(k << 1) + 1];
	z1 = c__[(i__ << 1) + 2] - c__[(k << 1) + 2];
	p[0] = x + ch[i__ + ch_dim1] * y - ch[i__ + (ch_dim1 << 1)] * z__;
	p[1] = z1 - ch[i__ + (ch_dim1 << 1)] * y - ch[i__ + ch_dim1] * z__;
	i__2 = i__;
	q__1.r = p[0], q__1.i = p[1];
	cy[i__2].r = q__1.r, cy[i__2].i = q__1.i;
	--k;
/* L10: */
    }
    return 0;
} /* rconv2_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     CFFT2 */
/*                                                                  ************ */
/* Subroutine */ int cfft2_(integer *init, integer *ix, integer *n, complex *
	cx, complex *cwork, complex *cy)
{
    static integer l2, n2, ls, ns;
    extern /* Subroutine */ int abort_(char *, ftnlen), cable2_(integer *, 
	    complex *), ctock1_(integer *, complex *, complex *), ctock2_(
	    integer *, integer *, complex *, complex *, complex *), ctock3_(
	    integer *, integer *, complex *, complex *, complex *);
    extern integer popcnt_(integer *);


/*     THE STOCKHAM AUTO-SORT FFT */


/*     IS N IS THE POWER OF 2 ? */

    /* Parameter adjustments */
    --cy;
    --cwork;
    --cx;

    /* Function Body */
    l2 = popcnt_(n);
    if (l2 != 1) {
	abort_("CFFT2 ", (ftnlen)6);
    }
    n2 = *n + *n + 1;

/*     WANT TO GENERATE THE TABLE FOR SINE AND COSINE, IF INIT .NE. 0. */

    ns = *n / 2;
    if (*init == 0) {
	goto L10;
    }

/*     IF N = 0 OR 2, RETURN. */
/*     IF N = 4,8 OR 16 THEN USING TABLE TO GENERATE.(CBALE1) */
/*     IF N > 16 THEN USING THE NEW METHOD TO GENERATE.(CBLE2) */

    if (*n < 4) {
	abort_("CFFT2 ", (ftnlen)6);
    }
    cable2_(&ns, &cwork[n2]);
    return 0;

/*     USING THE PREVIOUS TABLE TO PERFORM FOURIER TRANSFORMATION, */
/*     IF INIT = 0. */

L10:
    ctock1_(&ns, &cwork[1], &cx[1]);
    if (*ix <= 0) {
	goto L50;
    }

/*     THE FOLLOWING STATEMENTS FOR FOURIER ANALYSIS. */
/*            ^I.^E. IX > 0. */

    ls = 2;
    ns /= 2;
L20:
    if (ns == 1) {
	goto L120;
    }
    ctock2_(&ls, &ns, &cwork[*n + 1], &cwork[1], &cwork[n2]);
    ls += ls;
    ns /= 2;
    if (ns == 1) {
	goto L30;
    }
    ctock2_(&ls, &ns, &cwork[1], &cwork[*n + 1], &cwork[n2]);
    ls += ls;
    ns /= 2;
    goto L20;
L30:
    ctock2_(&ls, &ns, &cy[1], &cwork[*n + 1], &cwork[n2]);
    return 0;
L120:
    ctock2_(&ls, &ns, &cy[1], &cwork[1], &cwork[n2]);
    return 0;

/*     THE FOLLOWING STATEMENTS FOR FOURIER SYNTHESIS. */
/*            ^I.^E. IX < 0. */

L50:
    ls = 2;
    ns /= 2;
L60:
    if (ns == 1) {
	goto L130;
    }
    ctock3_(&ls, &ns, &cwork[*n + 1], &cwork[1], &cwork[n2]);
    ls += ls;
    ns /= 2;
    if (ns == 1) {
	goto L70;
    }
    ctock3_(&ls, &ns, &cwork[1], &cwork[*n + 1], &cwork[n2]);
    ls += ls;
    ns /= 2;
    goto L60;
L70:
    ctock3_(&ls, &ns, &cy[1], &cwork[*n + 1], &cwork[n2]);
    return 0;
L130:
    ctock3_(&ls, &ns, &cy[1], &cwork[1], &cwork[n2]);
    return 0;
} /* cfft2_ */

/* -----------------------------------------------  ************ */
/*                                                    CTOCK1 */
/*                                                 ************ */
/* Subroutine */ int ctock1_(integer *ns, complex *c__, complex *ch)
{
    /* System generated locals */
    integer c_dim1, c_offset, ch_dim1, ch_offset, i__1, i__2, i__3, i__4;
    complex q__1;

    /* Local variables */
    static integer j;



    /* Parameter adjustments */
    ch_dim1 = *ns;
    ch_offset = 1 + ch_dim1 * 3;
    ch -= ch_offset;
    c_dim1 = *ns;
    c_offset = 1 + (c_dim1 << 1);
    c__ -= c_offset;

    /* Function Body */
    i__1 = *ns;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + (c_dim1 << 1);
	i__3 = j + ch_dim1 * 3;
	i__4 = j + (ch_dim1 << 2);
	q__1.r = ch[i__3].r + ch[i__4].r, q__1.i = ch[i__3].i + ch[i__4].i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
	i__2 = j + c_dim1 * 3;
	i__3 = j + ch_dim1 * 3;
	i__4 = j + (ch_dim1 << 2);
	q__1.r = ch[i__3].r - ch[i__4].r, q__1.i = ch[i__3].i - ch[i__4].i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
/* L300: */
    }
    return 0;
} /* ctock1_ */

/* -----------------------------------------------  ************ */
/*                                                    CTOCK2 */
/*                                                 ************ */
/* Subroutine */ int ctock2_(integer *ls, integer *ns, complex *c__, complex *
	ch, complex *ch2)
{
    /* System generated locals */
    integer c_dim1, c_dim2, c_offset, ch_dim1, ch_offset, ch2_dim1, 
	    ch2_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    complex q__1, q__2;

    /* Local variables */
    static integer i__, j;



    /* Parameter adjustments */
    ch2_dim1 = *ns;
    ch2_offset = 1 + ch2_dim1;
    ch2 -= ch2_offset;
    ch_dim1 = *ns;
    ch_offset = 1 + ch_dim1 * 3;
    ch -= ch_offset;
    c_dim1 = *ns;
    c_dim2 = *ls;
    c_offset = 1 + c_dim1 * (1 + c_dim2);
    c__ -= c_offset;

    /* Function Body */
    if (*ls > *ns) {
	goto L20;
    }
    i__1 = *ls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ns;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    i__5 = i__ * ch2_dim1 + 1;
	    i__6 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__2.r = ch2[i__5].r * ch[i__6].r - ch2[i__5].i * ch[i__6].i, 
		    q__2.i = ch2[i__5].r * ch[i__6].i + ch2[i__5].i * ch[i__6]
		    .r;
	    q__1.r = ch[i__4].r + q__2.r, q__1.i = ch[i__4].i + q__2.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    i__5 = i__ * ch2_dim1 + 1;
	    i__6 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__2.r = ch2[i__5].r * ch[i__6].r - ch2[i__5].i * ch[i__6].i, 
		    q__2.i = ch2[i__5].r * ch[i__6].i + ch2[i__5].i * ch[i__6]
		    .r;
	    q__1.r = ch[i__4].r - q__2.r, q__1.i = ch[i__4].i - q__2.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L200: */
	}
    }
    return 0;
L20:
    i__2 = *ns;
    for (j = 1; j <= i__2; ++j) {
	i__1 = *ls;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    i__5 = i__ * ch2_dim1 + 1;
	    i__6 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__2.r = ch2[i__5].r * ch[i__6].r - ch2[i__5].i * ch[i__6].i, 
		    q__2.i = ch2[i__5].r * ch[i__6].i + ch2[i__5].i * ch[i__6]
		    .r;
	    q__1.r = ch[i__4].r + q__2.r, q__1.i = ch[i__4].i + q__2.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    i__5 = i__ * ch2_dim1 + 1;
	    i__6 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__2.r = ch2[i__5].r * ch[i__6].r - ch2[i__5].i * ch[i__6].i, 
		    q__2.i = ch2[i__5].r * ch[i__6].i + ch2[i__5].i * ch[i__6]
		    .r;
	    q__1.r = ch[i__4].r - q__2.r, q__1.i = ch[i__4].i - q__2.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L400: */
	}
    }
    return 0;
} /* ctock2_ */

/* -----------------------------------------------  ************ */
/*                                                    CTOCK3 */
/*                                                 ************ */
/* Subroutine */ int ctock3_(integer *ls, integer *ns, complex *c__, complex *
	ch, complex *ch2)
{
    /* System generated locals */
    integer c_dim1, c_dim2, c_offset, ch_dim1, ch_offset, ch2_dim1, 
	    ch2_offset, i__1, i__2, i__3, i__4, i__5;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    static integer i__, j;



    /* Parameter adjustments */
    ch2_dim1 = *ns;
    ch2_offset = 1 + ch2_dim1;
    ch2 -= ch2_offset;
    ch_dim1 = *ns;
    ch_offset = 1 + ch_dim1 * 3;
    ch -= ch_offset;
    c_dim1 = *ns;
    c_dim2 = *ls;
    c_offset = 1 + c_dim1 * (1 + c_dim2);
    c__ -= c_offset;

    /* Function Body */
    if (*ls > *ns) {
	goto L30;
    }
    i__1 = *ls;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ns;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    r_cnjg(&q__3, &ch2[i__ * ch2_dim1 + 1]);
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__2.r = q__3.r * ch[i__5].r - q__3.i * ch[i__5].i, q__2.i = 
		    q__3.r * ch[i__5].i + q__3.i * ch[i__5].r;
	    q__1.r = ch[i__4].r + q__2.r, q__1.i = ch[i__4].i + q__2.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    r_cnjg(&q__3, &ch2[i__ * ch2_dim1 + 1]);
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__2.r = q__3.r * ch[i__5].r - q__3.i * ch[i__5].i, q__2.i = 
		    q__3.r * ch[i__5].i + q__3.i * ch[i__5].r;
	    q__1.r = ch[i__4].r - q__2.r, q__1.i = ch[i__4].i - q__2.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L600: */
	}
    }
    return 0;
L30:
    i__2 = *ns;
    for (j = 1; j <= i__2; ++j) {
	i__1 = *ls;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__3 = j + (i__ + c_dim2) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    r_cnjg(&q__3, &ch2[i__ * ch2_dim1 + 1]);
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__2.r = q__3.r * ch[i__5].r - q__3.i * ch[i__5].i, q__2.i = 
		    q__3.r * ch[i__5].i + q__3.i * ch[i__5].r;
	    q__1.r = ch[i__4].r + q__2.r, q__1.i = ch[i__4].i + q__2.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
	    i__3 = j + (i__ + (c_dim2 << 1)) * c_dim1;
	    i__4 = j + ((i__ << 1) + 1) * ch_dim1;
	    r_cnjg(&q__3, &ch2[i__ * ch2_dim1 + 1]);
	    i__5 = j + ((i__ << 1) + 2) * ch_dim1;
	    q__2.r = q__3.r * ch[i__5].r - q__3.i * ch[i__5].i, q__2.i = 
		    q__3.r * ch[i__5].i + q__3.i * ch[i__5].r;
	    q__1.r = ch[i__4].r - q__2.r, q__1.i = ch[i__4].i - q__2.i;
	    c__[i__3].r = q__1.r, c__[i__3].i = q__1.i;
/* L800: */
	}
    }
    return 0;
} /* ctock3_ */

/* -----------------------------------------------  ************ */
/*                                                    CABLE2 */
/*                                                 ************ */
/* Subroutine */ int cable2_(integer *nn, real *work)
{
    /* Initialized data */

    static real twopi = 6.28318530717958647692f;

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, n;
    static real p2;


    /* Parameter adjustments */
    work -= 3;

    /* Function Body */

    n = 2 * *nn;
    p2 = twopi / n;
    i__1 = *nn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	work[(i__ << 1) + 1] = cos(p2 * (i__ - 1));
	work[(i__ << 1) + 2] = sin(p2 * (i__ - 1));
/* L10: */
    }
    return 0;
} /* cable2_ */

/* -------------------------------------------------------------     ************ */
/*                                                                      ABORT */
/*                                                                  ************ */
/* Subroutine */ int abort_(char *nme, ftnlen nme_len)
{
    /* Builtin functions */
    integer s_wsle(cilist *), e_wsle(void), do_lio(integer *, integer *, char 
	    *, ftnlen);

    /* Local variables */
    static real den, num, div0;

    /* Fortran I/O blocks */
    static cilist io___66 = { 0, 6, 0, 0, 0 };
    static cilist io___67 = { 0, 6, 0, 0, 0 };
    static cilist io___68 = { 0, 6, 0, 0, 0 };



/*     Routine to abort execution if N is not of the form 2**I. */
/*     Martin J. McBride.  2/27/86. */
/*     General Electric CRD, Information System Operation. */

    s_wsle(&io___66);
    e_wsle();
    s_wsle(&io___67);
    do_lio(&c__9, &c__1, nme, (ftnlen)6);
    do_lio(&c__9, &c__1, " called with N not of the form N=2**I where I=>2.", 
	    (ftnlen)49);
    e_wsle();
    s_wsle(&io___68);
    e_wsle();
    num = 1.f;
    den = 0.f;
    div0 = num / den;
    return 0;
} /* abort_ */

/* -------------------------------------------------------------     ************ */
/*                                                                     POPCNT */
/*                                                                  ************ */
integer popcnt_(integer *n)
{
    /* System generated locals */
    integer ret_val;

    /* Local variables */
    static integer i__, ib, tmp;


/*     Routine to determine how many bits of N are 1 (is N = 2**I). */
/*     Martin J. McBride.  3/4/86. */
/*     General Electric CRD, Information System Operation. */

    i__ = 0;
    tmp = *n;
L10:
    ib = tmp % 2;
    if (ib != 0) {
	++i__;
    }
    tmp /= 2;
    if (tmp != 0) {
	goto L10;
    }
    ret_val = i__;
    return ret_val;
} /* popcnt_ */

