/* linpackd.f -- translated by f2c (version 20100827).
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

static integer c__1 = 1;
static doublereal c_b1033 = -1.;
static doublecomplex c_b1092 = {1.,0.};
static doublecomplex c_b1136 = {-1.,-0.};
static doublecomplex c_b2561 = {-1.,0.};

/* Subroutine */ int sgeco_(doublereal *a, integer *lda, integer *n, integer *
	ipvt, doublereal *rcond, doublereal *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer j, k, l;
    static doublereal s, t;
    static integer kb;
    static doublereal ek, sm, wk;
    static integer kp1;
    static doublereal wkm;
    static integer info;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int sgefa_(doublereal *, integer *, integer *, 
	    integer *, integer *), sscal_(integer *, doublereal *, doublereal 
	    *, integer *);
    static doublereal anorm;
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DGECO FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION */
/*     AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW DGECO BY DGESL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW DGECO BY DGESL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW DGECO BY DGEDI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW DGECO BY DGEDI. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE MATRIX TO BE FACTORED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS */
/*                WHICH WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE */
/*                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER */
/*                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR. */

/*        IPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        RCOND   DOUBLE PRECISION */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       DOUBLE PRECISION(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK SGEFA */
/*     BLAS SAXPY,SDOT,SSCAL,SASUM */
/*     FORTRAN DABS,DMAX1,DSIGN */

/*     INTERNAL VARIABLES */



/*     COMPUTE 1-NORM OF A */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --z__;

    /* Function Body */
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = sasum_(n, &a[j * a_dim1 + 1], &c__1);
	anorm = max(d__1,d__2);
/* L10: */
    }

/*     FACTOR */

    sgefa_(&a[a_offset], lda, n, &ipvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E . */
/*     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE */
/*     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE */
/*     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID */
/*     OVERFLOW. */

/*     SOLVE TRANS(U)*W = E */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L20: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.) {
	    d__1 = -z__[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z__[k], abs(d__1)) <= (d__2 = a[k + k * a_dim1], abs(
		d__2))) {
	    goto L30;
	}
	s = (d__1 = a[k + k * a_dim1], abs(d__1)) / (d__2 = ek - z__[k], abs(
		d__2));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = abs(wk);
	sm = abs(wkm);
	if (a[k + k * a_dim1] == 0.) {
	    goto L40;
	}
	wk /= a[k + k * a_dim1];
	wkm /= a[k + k * a_dim1];
	goto L50;
L40:
	wk = 1.;
	wkm = 1.;
L50:
	kp1 = k + 1;
	if (kp1 > *n) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    sm += (d__1 = z__[j] + wkm * a[k + j * a_dim1], abs(d__1));
	    z__[j] += wk * a[k + j * a_dim1];
	    s += (d__1 = z__[j], abs(d__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	t = wkm - wk;
	wk = wkm;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    z__[j] += t * a[k + j * a_dim1];
/* L70: */
	}
L80:
L90:
	z__[k] = wk;
/* L100: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*     SOLVE TRANS(L)*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if (k < *n) {
	    i__2 = *n - k;
	    z__[k] += sdot_(&i__2, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1],
		     &c__1);
	}
	if ((d__1 = z__[k], abs(d__1)) <= 1.) {
	    goto L110;
	}
	s = 1. / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* L120: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE L*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
	if (k < *n) {
	    i__2 = *n - k;
	    saxpy_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1], &
		    c__1);
	}
	if ((d__1 = z__[k], abs(d__1)) <= 1.) {
	    goto L130;
	}
	s = 1. / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE  U*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= (d__2 = a[k + k * a_dim1], abs(d__2)
		)) {
	    goto L150;
	}
	s = (d__1 = a[k + k * a_dim1], abs(d__1)) / (d__2 = z__[k], abs(d__2))
		;
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	if (a[k + k * a_dim1] != 0.) {
	    z__[k] /= a[k + k * a_dim1];
	}
	if (a[k + k * a_dim1] == 0.) {
	    z__[k] = 1.;
	}
	t = -z__[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
/* L160: */
    }
/*     MAKE ZNORM = 1.0 */
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* sgeco_ */

/* Subroutine */ int sgefa_(doublereal *a, integer *lda, integer *n, integer *
	ipvt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, k, l;
    static doublereal t;
    static integer kp1, nm1;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), saxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);
    extern integer isamax_(integer *, doublereal *, integer *);


/*     SGEFA FACTORS A DOUBLE PRECISION MATRIX BY GAUSSIAN ELIMINATION. */

/*     SGEFA IS USUALLY CALLED BY DGECO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */
/*     (TIME FOR DGECO) = (1 + 9/N)*(TIME FOR SGEFA) . */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE MATRIX TO BE FACTORED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS */
/*                WHICH WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE */
/*                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER */
/*                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR. */

/*        IPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR */
/*                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES */
/*                     INDICATE THAT DGESL OR DGEDI WILL DIVIDE BY ZERO */
/*                     IF CALLED.  USE  RCOND  IN DGECO FOR A RELIABLE */
/*                     INDICATION OF SINGULARITY. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SSCAL,ISAMAX */

/*     INTERNAL VARIABLES */



/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;

    /* Function Body */
    *info = 0;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        FIND L = PIVOT INDEX */

	i__2 = *n - k + 1;
	l = isamax_(&i__2, &a[k + k * a_dim1], &c__1) + k - 1;
	ipvt[k] = l;

/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */

	if (a[l + k * a_dim1] == 0.) {
	    goto L40;
	}

/*           INTERCHANGE IF NECESSARY */

	if (l == k) {
	    goto L10;
	}
	t = a[l + k * a_dim1];
	a[l + k * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = t;
L10:

/*           COMPUTE MULTIPLIERS */

	t = -1. / a[k + k * a_dim1];
	i__2 = *n - k;
	sscal_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1);

/*           ROW ELIMINATION WITH COLUMN INDEXING */

	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[l + j * a_dim1];
	    if (l == k) {
		goto L20;
	    }
	    a[l + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = t;
L20:
	    i__3 = *n - k;
	    saxpy_(&i__3, &t, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 + j * 
		    a_dim1], &c__1);
/* L30: */
	}
	goto L50;
L40:
	*info = k;
L50:
/* L60: */
	;
    }
L70:
    ipvt[*n] = *n;
    if (a[*n + *n * a_dim1] == 0.) {
	*info = *n;
    }
    return 0;
} /* sgefa_ */

/* Subroutine */ int sgesl_(doublereal *a, integer *lda, integer *n, integer *
	ipvt, doublereal *b, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer k, l;
    static doublereal t;
    static integer kb, nm1;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DGESL SOLVES THE DOUBLE PRECISION SYSTEM */
/*     A * X = B  OR  TRANS(A) * X = B */
/*     USING THE FACTORS COMPUTED BY DGECO OR SGEFA. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE OUTPUT FROM DGECO OR SGEFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        IPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM DGECO OR SGEFA. */

/*        B       DOUBLE PRECISION(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*        JOB     INTEGER */
/*                = 0         TO SOLVE  A*X = B , */
/*                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE */
/*                            TRANS(A)  IS THE TRANSPOSE. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A */
/*        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY */
/*        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER */
/*        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE */
/*        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0 */
/*        OR SGEFA HAS SET INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL DGECO(A,LDA,N,IPVT,RCOND,Z) */
/*           IF (RCOND IS TOO SMALL) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DGESL(A,LDA,N,IPVT,C(1,J),0) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SDOT */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --b;

    /* Function Body */
    nm1 = *n - 1;
    if (*job != 0) {
	goto L50;
    }

/*        JOB = 0 , SOLVE  A * X = B */
/*        FIRST SOLVE  L*Y = B */

    if (nm1 < 1) {
	goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = b[l];
	if (l == k) {
	    goto L10;
	}
	b[l] = b[k];
	b[k] = t;
L10:
	i__2 = *n - k;
	saxpy_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &b[k + 1], &c__1);
/* L20: */
    }
L30:

/*        NOW SOLVE  U*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L40: */
    }
    goto L100;
L50:

/*        JOB = NONZERO, SOLVE  TRANS(A) * X = B */
/*        FIRST SOLVE  TRANS(U)*Y = B */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	t = sdot_(&i__2, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
	b[k] = (b[k] - t) / a[k + k * a_dim1];
/* L60: */
    }

/*        NOW SOLVE TRANS(L)*X = Y */

    if (nm1 < 1) {
	goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
	i__2 = *n - k;
	b[k] += sdot_(&i__2, &a[k + 1 + k * a_dim1], &c__1, &b[k + 1], &c__1);
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	t = b[l];
	b[l] = b[k];
	b[k] = t;
L70:
/* L80: */
	;
    }
L90:
L100:
    return 0;
} /* sgesl_ */

/* Subroutine */ int sgedi_(doublereal *a, integer *lda, integer *n, integer *
	ipvt, doublereal *det, doublereal *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal t;
    static integer kb, kp1, nm1;
    static doublereal ten;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sswap_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX */
/*     USING THE FACTORS COMPUTED BY DGECO OR SGEFA. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE OUTPUT FROM DGECO OR SGEFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        IPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM DGECO OR SGEFA. */

/*        WORK    DOUBLE PRECISION(N) */
/*                WORK VECTOR.  CONTENTS DESTROYED. */

/*        JOB     INTEGER */
/*                = 11   BOTH DETERMINANT AND INVERSE. */
/*                = 01   INVERSE ONLY. */
/*                = 10   DETERMINANT ONLY. */

/*     ON RETURN */

/*        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED. */
/*                OTHERWISE UNCHANGED. */

/*        DET     DOUBLE PRECISION(2) */
/*                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED. */
/*                OTHERWISE NOT REFERENCED. */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0 */
/*                OR  DET(1) .EQ. 0.0 . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED. */
/*        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY */
/*        AND IF DGECO HAS SET RCOND .GT. 0.0 OR SGEFA HAS SET */
/*        INFO .EQ. 0 . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SSCAL,SSWAP */
/*     FORTRAN DABS,MOD */

/*     INTERNAL VARIABLES */



/*     COMPUTE DETERMINANT */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --det;
    --work;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    det[1] = -det[1];
	}
	det[1] = a[i__ + i__ * a_dim1] * det[1];
/*        ...EXIT */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (abs(det[1]) >= 1.) {
	    goto L20;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (abs(det[1]) < ten) {
	    goto L40;
	}
	det[1] /= ten;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*     COMPUTE INVERSE(U) */

    if (*job % 10 == 0) {
	goto L150;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	a[k + k * a_dim1] = 1. / a[k + k * a_dim1];
	t = -a[k + k * a_dim1];
	i__2 = k - 1;
	sscal_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[k + j * a_dim1];
	    a[k + j * a_dim1] = 0.;
	    saxpy_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        FORM INVERSE(U)*INVERSE(L) */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L140;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
	kp1 = k + 1;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    work[i__] = a[i__ + k * a_dim1];
	    a[i__ + k * a_dim1] = 0.;
/* L110: */
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = work[j];
	    saxpy_(n, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L120: */
	}
	l = ipvt[k];
	if (l != k) {
	    sswap_(n, &a[k * a_dim1 + 1], &c__1, &a[l * a_dim1 + 1], &c__1);
	}
/* L130: */
    }
L140:
L150:
    return 0;
} /* sgedi_ */

/* Subroutine */ int sgbco_(doublereal *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublereal *rcond, 
	doublereal *z__)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer j, k, l, m;
    static doublereal s, t;
    static integer kb, la;
    static doublereal ek;
    static integer lm, mm, is, ju;
    static doublereal sm, wk;
    static integer lz, kp1;
    static doublereal wkm;
    static integer info;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int sgbfa_(doublereal *, integer *, integer *, 
	    integer *, integer *, integer *, integer *), sscal_(integer *, 
	    doublereal *, doublereal *, integer *);
    static doublereal anorm;
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DGBCO FACTORS A DOUBLE PRECISION BAND MATRIX BY GAUSSIAN */
/*     ELIMINATION AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, SGBFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW DGBCO BY DGBSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW DGBCO BY DGBSL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW DGBCO BY DGBDI. */

/*     ON ENTRY */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS */
/*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND */
/*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
/*                ML+1 THROUGH 2*ML+MU+1 OF  ABD . */
/*                SEE THE COMMENTS BELOW FOR DETAILS. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */
/*                LDA MUST BE .GE. 2*ML + MU + 1 . */

/*        N       INTEGER */
/*                THE ORDER OF THE ORIGINAL MATRIX. */

/*        ML      INTEGER */
/*                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL. */
/*                0 .LE. ML .LT. N . */

/*        MU      INTEGER */
/*                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */
/*                0 .LE. MU .LT. N . */
/*                MORE EFFICIENT IF  ML .LE. MU . */

/*     ON RETURN */

/*        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
/*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE */
/*                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER */
/*                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR. */

/*        IPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        RCOND   DOUBLE PRECISION */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       DOUBLE PRECISION(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     BAND STORAGE */

/*           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT */
/*           WILL SET UP THE INPUT. */

/*                   ML = (BAND WIDTH BELOW THE DIAGONAL) */
/*                   MU = (BAND WIDTH ABOVE THE DIAGONAL) */
/*                   M = ML + MU + 1 */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX0(1, J-MU) */
/*                      I2 = MIN0(N, J+ML) */
/*                      DO 10 I = I1, I2 */
/*                         K = I - J + M */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD . */
/*           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR */
/*           ELEMENTS GENERATED DURING THE TRIANGULARIZATION. */
/*           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 . */
/*           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE */
/*           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED. */

/*     EXAMPLE..  IF THE ORIGINAL MATRIX IS */

/*           11 12 13  0  0  0 */
/*           21 22 23 24  0  0 */
/*            0 32 33 34 35  0 */
/*            0  0 43 44 45 46 */
/*            0  0  0 54 55 56 */
/*            0  0  0  0 65 66 */

/*      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN */

/*            *  *  *  +  +  +  , * = NOT USED */
/*            *  * 13 24 35 46  , + = USED FOR PIVOTING */
/*            * 12 23 34 45 56 */
/*           11 22 33 44 55 66 */
/*           21 32 43 54 65  * */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK SGBFA */
/*     BLAS SAXPY,SDOT,SSCAL,SASUM */
/*     FORTRAN DABS,DMAX1,MAX0,MIN0,DSIGN */

/*     INTERNAL VARIABLES */



/*     COMPUTE 1-NORM OF A */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --z__;

    /* Function Body */
    anorm = 0.;
    l = *ml + 1;
    is = l + *mu;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = sasum_(&l, &abd[is + j * abd_dim1], &c__1);
	anorm = max(d__1,d__2);
	if (is > *ml + 1) {
	    --is;
	}
	if (j <= *mu) {
	    ++l;
	}
	if (j >= *n - *ml) {
	    --l;
	}
/* L10: */
    }

/*     FACTOR */

    sgbfa_(&abd[abd_offset], lda, n, ml, mu, &ipvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E . */
/*     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE */
/*     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE */
/*     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID */
/*     OVERFLOW. */

/*     SOLVE TRANS(U)*W = E */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L20: */
    }
    m = *ml + *mu + 1;
    ju = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.) {
	    d__1 = -z__[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z__[k], abs(d__1)) <= (d__2 = abd[m + k * abd_dim1], 
		abs(d__2))) {
	    goto L30;
	}
	s = (d__1 = abd[m + k * abd_dim1], abs(d__1)) / (d__2 = ek - z__[k], 
		abs(d__2));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = abs(wk);
	sm = abs(wkm);
	if (abd[m + k * abd_dim1] == 0.) {
	    goto L40;
	}
	wk /= abd[m + k * abd_dim1];
	wkm /= abd[m + k * abd_dim1];
	goto L50;
L40:
	wk = 1.;
	wkm = 1.;
L50:
	kp1 = k + 1;
/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = max(i__3,i__4);
	ju = min(i__2,*n);
	mm = m;
	if (kp1 > ju) {
	    goto L90;
	}
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --mm;
	    sm += (d__1 = z__[j] + wkm * abd[mm + j * abd_dim1], abs(d__1));
	    z__[j] += wk * abd[mm + j * abd_dim1];
	    s += (d__1 = z__[j], abs(d__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	t = wkm - wk;
	wk = wkm;
	mm = m;
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --mm;
	    z__[j] += t * abd[mm + j * abd_dim1];
/* L70: */
	}
L80:
L90:
	z__[k] = wk;
/* L100: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*     SOLVE TRANS(L)*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	if (k < *n) {
	    z__[k] += sdot_(&lm, &abd[m + 1 + k * abd_dim1], &c__1, &z__[k + 
		    1], &c__1);
	}
	if ((d__1 = z__[k], abs(d__1)) <= 1.) {
	    goto L110;
	}
	s = 1. / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* L120: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE L*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	if (k < *n) {
	    saxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &z__[k + 1], &
		    c__1);
	}
	if ((d__1 = z__[k], abs(d__1)) <= 1.) {
	    goto L130;
	}
	s = 1. / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE  U*Z = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= (d__2 = abd[m + k * abd_dim1], abs(
		d__2))) {
	    goto L150;
	}
	s = (d__1 = abd[m + k * abd_dim1], abs(d__1)) / (d__2 = z__[k], abs(
		d__2));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	if (abd[m + k * abd_dim1] != 0.) {
	    z__[k] /= abd[m + k * abd_dim1];
	}
	if (abd[m + k * abd_dim1] == 0.) {
	    z__[k] = 1.;
	}
	lm = min(k,m) - 1;
	la = m - lm;
	lz = k - lm;
	t = -z__[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lz], &c__1);
/* L160: */
    }
/*     MAKE ZNORM = 1.0 */
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* sgbco_ */

/* Subroutine */ int sgbfa_(doublereal *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal t;
    static integer i0, j0, j1, lm, mm, ju, jz, kp1, nm1;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), saxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);
    extern integer isamax_(integer *, doublereal *, integer *);


/*     SGBFA FACTORS A DOUBLE PRECISION BAND MATRIX BY ELIMINATION. */

/*     SGBFA IS USUALLY CALLED BY DGBCO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */

/*     ON ENTRY */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS */
/*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND */
/*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
/*                ML+1 THROUGH 2*ML+MU+1 OF  ABD . */
/*                SEE THE COMMENTS BELOW FOR DETAILS. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */
/*                LDA MUST BE .GE. 2*ML + MU + 1 . */

/*        N       INTEGER */
/*                THE ORDER OF THE ORIGINAL MATRIX. */

/*        ML      INTEGER */
/*                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL. */
/*                0 .LE. ML .LT. N . */

/*        MU      INTEGER */
/*                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */
/*                0 .LE. MU .LT. N . */
/*                MORE EFFICIENT IF  ML .LE. MU . */
/*     ON RETURN */

/*        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
/*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE */
/*                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER */
/*                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR. */

/*        IPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR */
/*                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES */
/*                     INDICATE THAT DGBSL WILL DIVIDE BY ZERO IF */
/*                     CALLED.  USE  RCOND  IN DGBCO FOR A RELIABLE */
/*                     INDICATION OF SINGULARITY. */

/*     BAND STORAGE */

/*           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT */
/*           WILL SET UP THE INPUT. */

/*                   ML = (BAND WIDTH BELOW THE DIAGONAL) */
/*                   MU = (BAND WIDTH ABOVE THE DIAGONAL) */
/*                   M = ML + MU + 1 */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX0(1, J-MU) */
/*                      I2 = MIN0(N, J+ML) */
/*                      DO 10 I = I1, I2 */
/*                         K = I - J + M */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD . */
/*           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR */
/*           ELEMENTS GENERATED DURING THE TRIANGULARIZATION. */
/*           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 . */
/*           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE */
/*           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SSCAL,ISAMAX */
/*     FORTRAN MAX0,MIN0 */

/*     INTERNAL VARIABLES */



    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;

    /* Function Body */
    m = *ml + *mu + 1;
    *info = 0;

/*     ZERO INITIAL FILL-IN COLUMNS */

    j0 = *mu + 2;
    j1 = min(*n,m) - 1;
    if (j1 < j0) {
	goto L30;
    }
    i__1 = j1;
    for (jz = j0; jz <= i__1; ++jz) {
	i0 = m + 1 - jz;
	i__2 = *ml;
	for (i__ = i0; i__ <= i__2; ++i__) {
	    abd[i__ + jz * abd_dim1] = 0.;
/* L10: */
	}
/* L20: */
    }
L30:
    jz = j1;
    ju = 0;

/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L130;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        ZERO NEXT FILL-IN COLUMN */

	++jz;
	if (jz > *n) {
	    goto L50;
	}
	if (*ml < 1) {
	    goto L50;
	}
	i__2 = *ml;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    abd[i__ + jz * abd_dim1] = 0.;
/* L40: */
	}
L50:

/*        FIND L = PIVOT INDEX */

/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	i__2 = lm + 1;
	l = isamax_(&i__2, &abd[m + k * abd_dim1], &c__1) + m - 1;
	ipvt[k] = l + k - m;

/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */

	if (abd[l + k * abd_dim1] == 0.) {
	    goto L100;
	}

/*           INTERCHANGE IF NECESSARY */

	if (l == m) {
	    goto L60;
	}
	t = abd[l + k * abd_dim1];
	abd[l + k * abd_dim1] = abd[m + k * abd_dim1];
	abd[m + k * abd_dim1] = t;
L60:

/*           COMPUTE MULTIPLIERS */

	t = -1. / abd[m + k * abd_dim1];
	sscal_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1);

/*           ROW ELIMINATION WITH COLUMN INDEXING */

/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = max(i__3,i__4);
	ju = min(i__2,*n);
	mm = m;
	if (ju < kp1) {
	    goto L90;
	}
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --l;
	    --mm;
	    t = abd[l + j * abd_dim1];
	    if (l == mm) {
		goto L70;
	    }
	    abd[l + j * abd_dim1] = abd[mm + j * abd_dim1];
	    abd[mm + j * abd_dim1] = t;
L70:
	    saxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &abd[mm + 1 + 
		    j * abd_dim1], &c__1);
/* L80: */
	}
L90:
	goto L110;
L100:
	*info = k;
L110:
/* L120: */
	;
    }
L130:
    ipvt[*n] = *n;
    if (abd[m + *n * abd_dim1] == 0.) {
	*info = *n;
    }
    return 0;
} /* sgbfa_ */

/* Subroutine */ int sgbsl_(doublereal *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublereal *b, integer *job)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer k, l, m;
    static doublereal t;
    static integer kb, la, lb, lm, nm1;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DGBSL SOLVES THE DOUBLE PRECISION BAND SYSTEM */
/*     A * X = B  OR  TRANS(A) * X = B */
/*     USING THE FACTORS COMPUTED BY DGBCO OR SGBFA. */

/*     ON ENTRY */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                THE OUTPUT FROM DGBCO OR SGBFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */

/*        N       INTEGER */
/*                THE ORDER OF THE ORIGINAL MATRIX. */

/*        ML      INTEGER */
/*                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL. */

/*        MU      INTEGER */
/*                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */

/*        IPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM DGBCO OR SGBFA. */

/*        B       DOUBLE PRECISION(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*        JOB     INTEGER */
/*                = 0         TO SOLVE  A*X = B , */
/*                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE */
/*                            TRANS(A)  IS THE TRANSPOSE. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A */
/*        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY */
/*        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER */
/*        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE */
/*        CALLED CORRECTLY AND IF DGBCO HAS SET RCOND .GT. 0.0 */
/*        OR SGBFA HAS SET INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z) */
/*           IF (RCOND IS TOO SMALL) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SDOT */
/*     FORTRAN MIN0 */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --b;

    /* Function Body */
    m = *mu + *ml + 1;
    nm1 = *n - 1;
    if (*job != 0) {
	goto L50;
    }

/*        JOB = 0 , SOLVE  A * X = B */
/*        FIRST SOLVE L*Y = B */

    if (*ml == 0) {
	goto L30;
    }
    if (nm1 < 1) {
	goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	l = ipvt[k];
	t = b[l];
	if (l == k) {
	    goto L10;
	}
	b[l] = b[k];
	b[k] = t;
L10:
	saxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &c__1);
/* L20: */
    }
L30:

/*        NOW SOLVE  U*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= abd[m + k * abd_dim1];
	lm = min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	t = -b[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/* L40: */
    }
    goto L100;
L50:

/*        JOB = NONZERO, SOLVE  TRANS(A) * X = B */
/*        FIRST SOLVE  TRANS(U)*Y = B */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	lm = min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	t = sdot_(&lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
	b[k] = (b[k] - t) / abd[m + k * abd_dim1];
/* L60: */
    }

/*        NOW SOLVE TRANS(L)*X = Y */

    if (*ml == 0) {
	goto L90;
    }
    if (nm1 < 1) {
	goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	b[k] += sdot_(&lm, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &
		c__1);
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	t = b[l];
	b[l] = b[k];
	b[k] = t;
L70:
/* L80: */
	;
    }
L90:
L100:
    return 0;
} /* sgbsl_ */

/* Subroutine */ int sgbdi_(doublereal *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublereal *det)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1;

    /* Local variables */
    static integer i__, m;
    static doublereal ten;


/*     DGBDI COMPUTES THE DETERMINANT OF A BAND MATRIX */
/*     USING THE FACTORS COMPUTED BY DGBCO OR SGBFA. */
/*     IF THE INVERSE IS NEEDED, USE DGBSL  N  TIMES. */

/*     ON ENTRY */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                THE OUTPUT FROM DGBCO OR SGBFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */

/*        N       INTEGER */
/*                THE ORDER OF THE ORIGINAL MATRIX. */

/*        ML      INTEGER */
/*                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL. */

/*        MU      INTEGER */
/*                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */

/*        IPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM DGBCO OR SGBFA. */

/*     ON RETURN */

/*        DET     DOUBLE PRECISION(2) */
/*                DETERMINANT OF ORIGINAL MATRIX. */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0 */
/*                OR  DET(1) = 0.0 . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     FORTRAN DABS */

/*     INTERNAL VARIABLES */



    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --det;

    /* Function Body */
    m = *ml + *mu + 1;
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    det[1] = -det[1];
	}
	det[1] = abd[m + i__ * abd_dim1] * det[1];
/*     ...EXIT */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (abs(det[1]) >= 1.) {
	    goto L20;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (abs(det[1]) < ten) {
	    goto L40;
	}
	det[1] /= ten;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
    return 0;
} /* sgbdi_ */

/* Subroutine */ int spoco_(doublereal *a, integer *lda, integer *n, 
	doublereal *rcond, doublereal *z__, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, t;
    static integer kb;
    static doublereal ek, sm, wk;
    static integer jm1, kp1;
    static doublereal wkm;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), spofa_(doublereal *, integer *, integer *, integer *);
    static doublereal anorm;
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DPOCO FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE */
/*     MATRIX AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, SPOFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW DPOCO BY DPOSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW DPOCO BY DPOSL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW DPOCO BY DPODI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW DPOCO BY DPODI. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE */
/*                DIAGONAL AND UPPER TRIANGLE ARE USED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R */
/*                WHERE  TRANS(R)  IS THE TRANSPOSE. */
/*                THE STRICT LOWER TRIANGLE IS UNALTERED. */
/*                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE. */

/*        RCOND   DOUBLE PRECISION */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED. */

/*        Z       DOUBLE PRECISION(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */
/*                IF  INFO .NE. 0 , Z  IS UNCHANGED. */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR */
/*                     OF ORDER  K  IS NOT POSITIVE DEFINITE. */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK SPOFA */
/*     BLAS SAXPY,SDOT,SSCAL,SASUM */
/*     FORTRAN DABS,DMAX1,DREAL,DSIGN */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A USING ONLY UPPER HALF */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = sasum_(&j, &a[j * a_dim1 + 1], &c__1);
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = z__[j];
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    spofa_(&a[a_offset], lda, n, info);
    if (*info != 0) {
	goto L180;
    }

/*        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E . */
/*        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*        SOLVE TRANS(R)*W = E */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L50: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.) {
	    d__1 = -z__[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z__[k], abs(d__1)) <= a[k + k * a_dim1]) {
	    goto L60;
	}
	s = a[k + k * a_dim1] / (d__1 = ek - z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L60:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = abs(wk);
	sm = abs(wkm);
	wk /= a[k + k * a_dim1];
	wkm /= a[k + k * a_dim1];
	kp1 = k + 1;
	if (kp1 > *n) {
	    goto L100;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    sm += (d__1 = z__[j] + wkm * a[k + j * a_dim1], abs(d__1));
	    z__[j] += wk * a[k + j * a_dim1];
	    s += (d__1 = z__[j], abs(d__1));
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	t = wkm - wk;
	wk = wkm;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    z__[j] += t * a[k + j * a_dim1];
/* L80: */
	}
L90:
L100:
	z__[k] = wk;
/* L110: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*        SOLVE R*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= a[k + k * a_dim1]) {
	    goto L120;
	}
	s = a[k + k * a_dim1] / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
L120:
	z__[k] /= a[k + k * a_dim1];
	t = -z__[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
/* L130: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*        SOLVE TRANS(R)*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	z__[k] -= sdot_(&i__2, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
	if ((d__1 = z__[k], abs(d__1)) <= a[k + k * a_dim1]) {
	    goto L140;
	}
	s = a[k + k * a_dim1] / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	z__[k] /= a[k + k * a_dim1];
/* L150: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        SOLVE R*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= a[k + k * a_dim1]) {
	    goto L160;
	}
	s = a[k + k * a_dim1] / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	z__[k] /= a[k + k * a_dim1];
	t = -z__[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
/* L170: */
    }
/*        MAKE ZNORM = 1.0 */
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
L180:
    return 0;
} /* spoco_ */

/* Subroutine */ int spofa_(doublereal *a, integer *lda, integer *n, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal s, t;
    static integer jm1;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);


/*     SPOFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE */
/*     MATRIX. */

/*     SPOFA IS USUALLY CALLED BY DPOCO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */
/*     (TIME FOR DPOCO) = (1 + 18/N)*(TIME FOR SPOFA) . */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE */
/*                DIAGONAL AND UPPER TRIANGLE ARE USED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R */
/*                WHERE  TRANS(R)  IS THE TRANSPOSE. */
/*                THE STRICT LOWER TRIANGLE IS UNALTERED. */
/*                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE. */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR */
/*                     OF ORDER  K  IS NOT POSITIVE DEFINITE. */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SDOT */
/*     FORTRAN DSQRT */

/*     INTERNAL VARIABLES */

/*     BEGIN BLOCK WITH ...EXITS TO 40 */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = k - 1;
	    t = a[k + j * a_dim1] - sdot_(&i__3, &a[k * a_dim1 + 1], &c__1, &
		    a[j * a_dim1 + 1], &c__1);
	    t /= a[k + k * a_dim1];
	    a[k + j * a_dim1] = t;
	    s += t * t;
/* L10: */
	}
L20:
	s = a[j + j * a_dim1] - s;
/*     ......EXIT */
	if (s <= 0.) {
	    goto L40;
	}
	a[j + j * a_dim1] = sqrt(s);
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* spofa_ */

/* Subroutine */ int sposl_(doublereal *a, integer *lda, integer *n, 
	doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static integer k;
    static doublereal t;
    static integer kb;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DPOSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE */
/*     SYSTEM A * X = B */
/*     USING THE FACTORS COMPUTED BY DPOCO OR SPOFA. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE OUTPUT FROM DPOCO OR SPOFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        B       DOUBLE PRECISION(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES */
/*        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE */
/*        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED */
/*        CORRECTLY AND  INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL DPOCO(A,LDA,N,RCOND,Z,INFO) */
/*           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DPOSL(A,LDA,N,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SDOT */

/*     INTERNAL VARIABLES */


/*     SOLVE TRANS(R)*Y = B */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	t = sdot_(&i__2, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
	b[k] = (b[k] - t) / a[k + k * a_dim1];
/* L10: */
    }

/*     SOLVE R*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L20: */
    }
    return 0;
} /* sposl_ */

/* Subroutine */ int spodi_(doublereal *a, integer *lda, integer *n, 
	doublereal *det, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, t;
    static integer jm1, kp1;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), saxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);


/*     DPODI COMPUTES THE DETERMINANT AND INVERSE OF A CERTAIN */
/*     DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE MATRIX (SEE BELOW) */
/*     USING THE FACTORS COMPUTED BY DPOCO, SPOFA OR DQRDC. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE OUTPUT  A  FROM DPOCO OR SPOFA */
/*                OR THE OUTPUT  X  FROM DQRDC. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        JOB     INTEGER */
/*                = 11   BOTH DETERMINANT AND INVERSE. */
/*                = 01   INVERSE ONLY. */
/*                = 10   DETERMINANT ONLY. */

/*     ON RETURN */

/*        A       IF DPOCO OR SPOFA WAS USED TO FACTOR  A  THEN */
/*                DPODI PRODUCES THE UPPER HALF OF INVERSE(A) . */
/*                IF DQRDC WAS USED TO DECOMPOSE  X  THEN */
/*                DPODI PRODUCES THE UPPER HALF OF INVERSE(TRANS(X)*X) */
/*                WHERE TRANS(X) IS THE TRANSPOSE. */
/*                ELEMENTS OF  A  BELOW THE DIAGONAL ARE UNCHANGED. */
/*                IF THE UNITS DIGIT OF JOB IS ZERO,  A  IS UNCHANGED. */

/*        DET     DOUBLE PRECISION(2) */
/*                DETERMINANT OF  A  OR OF  TRANS(X)*X  IF REQUESTED. */
/*                OTHERWISE NOT REFERENCED. */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. DET(1) .LT. 10.0 */
/*                OR  DET(1) .EQ. 0.0 . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED. */
/*        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY */
/*        AND IF DPOCO OR SPOFA HAS SET INFO .EQ. 0 . */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SSCAL */
/*     FORTRAN MOD */

/*     INTERNAL VARIABLES */


/*     COMPUTE DETERMINANT */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --det;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    s = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = a[i__ + i__ * a_dim1];
	det[1] = d__1 * d__1 * det[1];
/*        ...EXIT */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (det[1] >= 1.) {
	    goto L20;
	}
	det[1] = s * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (det[1] < s) {
	    goto L40;
	}
	det[1] /= s;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*     COMPUTE INVERSE(R) */

    if (*job % 10 == 0) {
	goto L140;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	a[k + k * a_dim1] = 1. / a[k + k * a_dim1];
	t = -a[k + k * a_dim1];
	i__2 = k - 1;
	sscal_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[k + j * a_dim1];
	    a[k + j * a_dim1] = 0.;
	    saxpy_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        FORM  INVERSE(R) * TRANS(INVERSE(R)) */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L120;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    t = a[k + j * a_dim1];
	    saxpy_(&k, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L110: */
	}
L120:
	t = a[j + j * a_dim1];
	sscal_(&j, &t, &a[j * a_dim1 + 1], &c__1);
/* L130: */
    }
L140:
    return 0;
} /* spodi_ */

/* Subroutine */ int sppco_(doublereal *ap, integer *n, doublereal *rcond, 
	doublereal *z__, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, t;
    static integer j1, kb;
    static doublereal ek;
    static integer ij, kj, kk;
    static doublereal sm, wk;
    static integer jm1, kp1;
    static doublereal wkm;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), sppfa_(doublereal *, integer *, integer *);
    static doublereal anorm;
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DPPCO FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE */
/*     MATRIX STORED IN PACKED FORM */
/*     AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, SPPFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW DPPCO BY DPPSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW DPPCO BY DPPSL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW DPPCO BY DPPDI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW DPPCO BY DPPDI. */

/*     ON ENTRY */

/*        AP      DOUBLE PRECISION (N*(N+1)/2) */
/*                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE */
/*                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY */
/*                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 . */
/*                SEE COMMENTS BELOW FOR DETAILS. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        AP      AN UPPER TRIANGULAR MATRIX  R , STORED IN PACKED */
/*                FORM, SO THAT  A = TRANS(R)*R . */
/*                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE. */

/*        RCOND   DOUBLE PRECISION */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED. */

/*        Z       DOUBLE PRECISION(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */
/*                IF  INFO .NE. 0 , Z  IS UNCHANGED. */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR */
/*                     OF ORDER  K  IS NOT POSITIVE DEFINITE. */

/*     PACKED STORAGE */

/*          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER */
/*          TRIANGLE OF A SYMMETRIC MATRIX. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K) = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK SPPFA */
/*     BLAS SAXPY,SDOT,SSCAL,SASUM */
/*     FORTRAN DABS,DMAX1,DREAL,DSIGN */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A */

    /* Parameter adjustments */
    --z__;
    --ap;

    /* Function Body */
    j1 = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = sasum_(&j, &ap[j1], &c__1);
	ij = j1;
	j1 += j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] += (d__1 = ap[ij], abs(d__1));
	    ++ij;
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = z__[j];
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    sppfa_(&ap[1], n, info);
    if (*info != 0) {
	goto L180;
    }

/*        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E . */
/*        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*        SOLVE TRANS(R)*W = E */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L50: */
    }
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk += k;
	if (z__[k] != 0.) {
	    d__1 = -z__[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z__[k], abs(d__1)) <= ap[kk]) {
	    goto L60;
	}
	s = ap[kk] / (d__1 = ek - z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L60:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = abs(wk);
	sm = abs(wkm);
	wk /= ap[kk];
	wkm /= ap[kk];
	kp1 = k + 1;
	kj = kk + k;
	if (kp1 > *n) {
	    goto L100;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    sm += (d__1 = z__[j] + wkm * ap[kj], abs(d__1));
	    z__[j] += wk * ap[kj];
	    s += (d__1 = z__[j], abs(d__1));
	    kj += j;
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	t = wkm - wk;
	wk = wkm;
	kj = kk + k;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    z__[j] += t * ap[kj];
	    kj += j;
/* L80: */
	}
L90:
L100:
	z__[k] = wk;
/* L110: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*        SOLVE R*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= ap[kk]) {
	    goto L120;
	}
	s = ap[kk] / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
L120:
	z__[k] /= ap[kk];
	kk -= k;
	t = -z__[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
/* L130: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*        SOLVE TRANS(R)*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	z__[k] -= sdot_(&i__2, &ap[kk + 1], &c__1, &z__[1], &c__1);
	kk += k;
	if ((d__1 = z__[k], abs(d__1)) <= ap[kk]) {
	    goto L140;
	}
	s = ap[kk] / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	z__[k] /= ap[kk];
/* L150: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        SOLVE R*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= ap[kk]) {
	    goto L160;
	}
	s = ap[kk] / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	z__[k] /= ap[kk];
	kk -= k;
	t = -z__[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
/* L170: */
    }
/*        MAKE ZNORM = 1.0 */
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
L180:
    return 0;
} /* sppco_ */

/* Subroutine */ int sppfa_(doublereal *ap, integer *n, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal s, t;
    static integer jj, kj, kk, jm1;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);


/*     SPPFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE */
/*     MATRIX STORED IN PACKED FORM. */

/*     SPPFA IS USUALLY CALLED BY DPPCO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */
/*     (TIME FOR DPPCO) = (1 + 18/N)*(TIME FOR SPPFA) . */

/*     ON ENTRY */

/*        AP      DOUBLE PRECISION (N*(N+1)/2) */
/*                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE */
/*                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY */
/*                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 . */
/*                SEE COMMENTS BELOW FOR DETAILS. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        AP      AN UPPER TRIANGULAR MATRIX  R , STORED IN PACKED */
/*                FORM, SO THAT  A = TRANS(R)*R . */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT */
/*                     POSITIVE DEFINITE. */


/*     PACKED STORAGE */

/*          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER */
/*          TRIANGLE OF A SYMMETRIC MATRIX. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K) = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SDOT */
/*     FORTRAN DSQRT */

/*     INTERNAL VARIABLES */

/*     BEGIN BLOCK WITH ...EXITS TO 40 */


    /* Parameter adjustments */
    --ap;

    /* Function Body */
    jj = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	kj = jj;
	kk = 0;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    ++kj;
	    i__3 = k - 1;
	    t = ap[kj] - sdot_(&i__3, &ap[kk + 1], &c__1, &ap[jj + 1], &c__1);
	    kk += k;
	    t /= ap[kk];
	    ap[kj] = t;
	    s += t * t;
/* L10: */
	}
L20:
	jj += j;
	s = ap[jj] - s;
/*     ......EXIT */
	if (s <= 0.) {
	    goto L40;
	}
	ap[jj] = sqrt(s);
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* sppfa_ */

/* Subroutine */ int sppsl_(doublereal *ap, integer *n, doublereal *b)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer k;
    static doublereal t;
    static integer kb, kk;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DPPSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE */
/*     SYSTEM A * X = B */
/*     USING THE FACTORS COMPUTED BY DPPCO OR SPPFA. */

/*     ON ENTRY */

/*        AP      DOUBLE PRECISION (N*(N+1)/2) */
/*                THE OUTPUT FROM DPPCO OR SPPFA. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        B       DOUBLE PRECISION(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES */
/*        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE */
/*        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED */
/*        CORRECTLY AND  INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL DPPCO(AP,N,RCOND,Z,INFO) */
/*           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DPPSL(AP,N,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SDOT */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    --b;
    --ap;

    /* Function Body */
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	t = sdot_(&i__2, &ap[kk + 1], &c__1, &b[1], &c__1);
	kk += k;
	b[k] = (b[k] - t) / ap[kk];
/* L10: */
    }
    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= ap[kk];
	kk -= k;
	t = -b[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &ap[kk + 1], &c__1, &b[1], &c__1);
/* L20: */
    }
    return 0;
} /* sppsl_ */

/* Subroutine */ int sppdi_(doublereal *ap, integer *n, doublereal *det, 
	integer *job)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, t;
    static integer j1, k1, ii, jj, kj, kk, jm1, kp1;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), saxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);


/*     DPPDI COMPUTES THE DETERMINANT AND INVERSE */
/*     OF A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE MATRIX */
/*     USING THE FACTORS COMPUTED BY DPPCO OR SPPFA . */

/*     ON ENTRY */

/*        AP      DOUBLE PRECISION (N*(N+1)/2) */
/*                THE OUTPUT FROM DPPCO OR SPPFA. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        JOB     INTEGER */
/*                = 11   BOTH DETERMINANT AND INVERSE. */
/*                = 01   INVERSE ONLY. */
/*                = 10   DETERMINANT ONLY. */

/*     ON RETURN */

/*        AP      THE UPPER TRIANGULAR HALF OF THE INVERSE . */
/*                THE STRICT LOWER TRIANGLE IS UNALTERED. */

/*        DET     DOUBLE PRECISION(2) */
/*                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED. */
/*                OTHERWISE NOT REFERENCED. */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. DET(1) .LT. 10.0 */
/*                OR  DET(1) .EQ. 0.0 . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED. */
/*        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY */
/*        AND IF DPOCO OR SPOFA HAS SET INFO .EQ. 0 . */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SSCAL */
/*     FORTRAN MOD */

/*     INTERNAL VARIABLES */


/*     COMPUTE DETERMINANT */

    /* Parameter adjustments */
    --det;
    --ap;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    s = 10.;
    ii = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii += i__;
/* Computing 2nd power */
	d__1 = ap[ii];
	det[1] = d__1 * d__1 * det[1];
/*        ...EXIT */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (det[1] >= 1.) {
	    goto L20;
	}
	det[1] = s * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (det[1] < s) {
	    goto L40;
	}
	det[1] /= s;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*     COMPUTE INVERSE(R) */

    if (*job % 10 == 0) {
	goto L140;
    }
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	k1 = kk + 1;
	kk += k;
	ap[kk] = 1. / ap[kk];
	t = -ap[kk];
	i__2 = k - 1;
	sscal_(&i__2, &t, &ap[k1], &c__1);
	kp1 = k + 1;
	j1 = kk + 1;
	kj = kk + k;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = ap[kj];
	    ap[kj] = 0.;
	    saxpy_(&k, &t, &ap[k1], &c__1, &ap[j1], &c__1);
	    j1 += j;
	    kj += j;
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        FORM  INVERSE(R) * TRANS(INVERSE(R)) */

    jj = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	j1 = jj + 1;
	jj += j;
	jm1 = j - 1;
	k1 = 1;
	kj = j1;
	if (jm1 < 1) {
	    goto L120;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    t = ap[kj];
	    saxpy_(&k, &t, &ap[j1], &c__1, &ap[k1], &c__1);
	    k1 += k;
	    ++kj;
/* L110: */
	}
L120:
	t = ap[jj];
	sscal_(&j, &t, &ap[j1], &c__1);
/* L130: */
    }
L140:
    return 0;
} /* sppdi_ */

/* Subroutine */ int spbco_(doublereal *abd, integer *lda, integer *n, 
	integer *m, doublereal *rcond, doublereal *z__, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal s, t;
    static integer j2, kb, la, lb;
    static doublereal ek;
    static integer lm;
    static doublereal sm, wk;
    static integer mu, kp1;
    static doublereal wkm;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int spbfa_(doublereal *, integer *, integer *, 
	    integer *, integer *), sscal_(integer *, doublereal *, doublereal 
	    *, integer *);
    static doublereal anorm;
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DPBCO FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE */
/*     MATRIX STORED IN BAND FORM AND ESTIMATES THE CONDITION OF THE */
/*     MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, SPBFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW DPBCO BY DPBSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW DPBCO BY DPBSL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW DPBCO BY DPBDI. */

/*     ON ENTRY */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER */
/*                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE */
/*                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE */
/*                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */
/*                LDA MUST BE .GE. M + 1 . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        M       INTEGER */
/*                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */
/*                0 .LE. M .LT. N . */

/*     ON RETURN */

/*        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND */
/*                FORM, SO THAT  A = TRANS(R)*R . */
/*                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE. */

/*        RCOND   DOUBLE PRECISION */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED. */

/*        Z       DOUBLE PRECISION(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */
/*                IF  INFO .NE. 0 , Z  IS UNCHANGED. */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR */
/*                     OF ORDER  K  IS NOT POSITIVE DEFINITE. */

/*     BAND STORAGE */

/*           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX, */
/*           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT. */

/*                   M = (BAND WIDTH ABOVE DIAGONAL) */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX0(1, J-M) */
/*                      DO 10 I = I1, J */
/*                         K = I-J+M+1 */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           THIS USES  M + 1  ROWS OF  A , EXCEPT FOR THE  M BY M */
/*           UPPER LEFT TRIANGLE, WHICH IS IGNORED. */

/*     EXAMPLE..  IF THE ORIGINAL MATRIX IS */

/*           11 12 13  0  0  0 */
/*           12 22 23 24  0  0 */
/*           13 23 33 34 35  0 */
/*            0 24 34 44 45 46 */
/*            0  0 35 45 55 56 */
/*            0  0  0 46 56 66 */

/*     THEN  N = 6 , M = 2  AND  ABD  SHOULD CONTAIN */

/*            *  * 13 24 35 46 */
/*            * 12 23 34 45 56 */
/*           11 22 33 44 55 66 */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK SPBFA */
/*     BLAS SAXPY,SDOT,SSCAL,SASUM */
/*     FORTRAN DABS,DMAX1,MAX0,MIN0,DREAL,DSIGN */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	i__2 = j, i__3 = *m + 1;
	l = min(i__2,i__3);
/* Computing MAX */
	i__2 = *m + 2 - j;
	mu = max(i__2,1);
	z__[j] = sasum_(&l, &abd[mu + j * abd_dim1], &c__1);
	k = j - l;
	if (*m < mu) {
	    goto L20;
	}
	i__2 = *m;
	for (i__ = mu; i__ <= i__2; ++i__) {
	    ++k;
	    z__[k] += (d__1 = abd[i__ + j * abd_dim1], abs(d__1));
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = z__[j];
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    spbfa_(&abd[abd_offset], lda, n, m, info);
    if (*info != 0) {
	goto L180;
    }

/*        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E . */
/*        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*        SOLVE TRANS(R)*W = E */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L50: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.) {
	    d__1 = -z__[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z__[k], abs(d__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L60;
	}
	s = abd[*m + 1 + k * abd_dim1] / (d__1 = ek - z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L60:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = abs(wk);
	sm = abs(wkm);
	wk /= abd[*m + 1 + k * abd_dim1];
	wkm /= abd[*m + 1 + k * abd_dim1];
	kp1 = k + 1;
/* Computing MIN */
	i__2 = k + *m;
	j2 = min(i__2,*n);
	i__ = *m + 1;
	if (kp1 > j2) {
	    goto L100;
	}
	i__2 = j2;
	for (j = kp1; j <= i__2; ++j) {
	    --i__;
	    sm += (d__1 = z__[j] + wkm * abd[i__ + j * abd_dim1], abs(d__1));
	    z__[j] += wk * abd[i__ + j * abd_dim1];
	    s += (d__1 = z__[j], abs(d__1));
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	t = wkm - wk;
	wk = wkm;
	i__ = *m + 1;
	i__2 = j2;
	for (j = kp1; j <= i__2; ++j) {
	    --i__;
	    z__[j] += t * abd[i__ + j * abd_dim1];
/* L80: */
	}
L90:
L100:
	z__[k] = wk;
/* L110: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*        SOLVE  R*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L120;
	}
	s = abd[*m + 1 + k * abd_dim1] / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
L120:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	t = -z__[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L130: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*        SOLVE TRANS(R)*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	z__[k] -= sdot_(&lm, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
	if ((d__1 = z__[k], abs(d__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L140;
	}
	s = abd[*m + 1 + k * abd_dim1] / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* L150: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        SOLVE  R*Z = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L160;
	}
	s = abd[*m + 1 + k * abd_dim1] / (d__1 = z__[k], abs(d__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	t = -z__[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L170: */
    }
/*        MAKE ZNORM = 1.0 */
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
L180:
    return 0;
} /* spbco_ */

/* Subroutine */ int spbfa_(doublereal *abd, integer *lda, integer *n, 
	integer *m, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal s, t;
    static integer ik, jk, mu;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);


/*     SPBFA FACTORS A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE */
/*     MATRIX STORED IN BAND FORM. */

/*     SPBFA IS USUALLY CALLED BY DPBCO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */

/*     ON ENTRY */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER */
/*                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE */
/*                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE */
/*                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */
/*                LDA MUST BE .GE. M + 1 . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        M       INTEGER */
/*                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */
/*                0 .LE. M .LT. N . */

/*     ON RETURN */

/*        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND */
/*                FORM, SO THAT  A = TRANS(R)*R . */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT */
/*                     POSITIVE DEFINITE. */

/*     BAND STORAGE */

/*           IF  A  IS A SYMMETRIC POSITIVE DEFINITE BAND MATRIX, */
/*           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT. */

/*                   M = (BAND WIDTH ABOVE DIAGONAL) */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX0(1, J-M) */
/*                      DO 10 I = I1, J */
/*                         K = I-J+M+1 */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SDOT */
/*     FORTRAN MAX0,DSQRT */

/*     INTERNAL VARIABLES */

/*     BEGIN BLOCK WITH ...EXITS TO 40 */


    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	ik = *m + 1;
/* Computing MAX */
	i__2 = j - *m;
	jk = max(i__2,1);
/* Computing MAX */
	i__2 = *m + 2 - j;
	mu = max(i__2,1);
	if (*m < mu) {
	    goto L20;
	}
	i__2 = *m;
	for (k = mu; k <= i__2; ++k) {
	    i__3 = k - mu;
	    t = abd[k + j * abd_dim1] - sdot_(&i__3, &abd[ik + jk * abd_dim1],
		     &c__1, &abd[mu + j * abd_dim1], &c__1);
	    t /= abd[*m + 1 + jk * abd_dim1];
	    abd[k + j * abd_dim1] = t;
	    s += t * t;
	    --ik;
	    ++jk;
/* L10: */
	}
L20:
	s = abd[*m + 1 + j * abd_dim1] - s;
/*     ......EXIT */
	if (s <= 0.) {
	    goto L40;
	}
	abd[*m + 1 + j * abd_dim1] = sqrt(s);
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* spbfa_ */

/* Subroutine */ int spbsl_(doublereal *abd, integer *lda, integer *n, 
	integer *m, doublereal *b)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2;

    /* Local variables */
    static integer k;
    static doublereal t;
    static integer kb, la, lb, lm;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DPBSL SOLVES THE DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE */
/*     BAND SYSTEM  A*X = B */
/*     USING THE FACTORS COMPUTED BY DPBCO OR SPBFA. */

/*     ON ENTRY */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                THE OUTPUT FROM DPBCO OR SPBFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        M       INTEGER */
/*                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */

/*        B       DOUBLE PRECISION(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES */
/*        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE */
/*        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED */
/*        CORRECTLY AND  INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL DPBCO(ABD,LDA,N,RCOND,Z,INFO) */
/*           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DPBSL(ABD,LDA,N,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SDOT */
/*     FORTRAN MIN0 */

/*     INTERNAL VARIABLES */


/*     SOLVE TRANS(R)*Y = B */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	t = sdot_(&lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
	b[k] = (b[k] - t) / abd[*m + 1 + k * abd_dim1];
/* L10: */
    }

/*     SOLVE R*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	b[k] /= abd[*m + 1 + k * abd_dim1];
	t = -b[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/* L20: */
    }
    return 0;
} /* spbsl_ */

/* Subroutine */ int spbdi_(doublereal *abd, integer *lda, integer *n, 
	integer *m, doublereal *det)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal s;


/*     DPBDI COMPUTES THE DETERMINANT */
/*     OF A DOUBLE PRECISION SYMMETRIC POSITIVE DEFINITE BAND MATRIX */
/*     USING THE FACTORS COMPUTED BY DPBCO OR SPBFA. */
/*     IF THE INVERSE IS NEEDED, USE DPBSL  N  TIMES. */

/*     ON ENTRY */

/*        ABD     DOUBLE PRECISION(LDA, N) */
/*                THE OUTPUT FROM DPBCO OR SPBFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        M       INTEGER */
/*                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */

/*     ON RETURN */

/*        DET     DOUBLE PRECISION(2) */
/*                DETERMINANT OF ORIGINAL MATRIX IN THE FORM */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. DET(1) .LT. 10.0 */
/*                OR  DET(1) .EQ. 0.0 . */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */


/*     INTERNAL VARIABLES */


/*     COMPUTE DETERMINANT */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --det;

    /* Function Body */
    det[1] = 1.;
    det[2] = 0.;
    s = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = abd[*m + 1 + i__ * abd_dim1];
	det[1] = d__1 * d__1 * det[1];
/*     ...EXIT */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (det[1] >= 1.) {
	    goto L20;
	}
	det[1] = s * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (det[1] < s) {
	    goto L40;
	}
	det[1] /= s;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
    return 0;
} /* spbdi_ */

/* Subroutine */ int ssico_(doublereal *a, integer *lda, integer *n, integer *
	kpvt, doublereal *rcond, doublereal *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, t, ak, bk, ek;
    static integer kp, ks, jm1, kps;
    static doublereal akm1, bkm1;
    static integer info;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal denom;
    extern /* Subroutine */ int ssifa_(doublereal *, integer *, integer *, 
	    integer *, integer *), sscal_(integer *, doublereal *, doublereal 
	    *, integer *);
    static doublereal anorm;
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DSICO FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX BY ELIMINATION */
/*     WITH SYMMETRIC PIVOTING AND ESTIMATES THE CONDITION OF THE */
/*     MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, SSIFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW DSICO BY DSISL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW DSICO BY DSISL. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW DSICO BY DSIDI. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW DSICO BY DSIDI. */
/*     TO COMPUTE  INERTIA(A), FOLLOW DSICO BY DSIDI. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                THE SYMMETRIC MATRIX TO BE FACTORED. */
/*                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     OUTPUT */

/*        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE */
/*                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        RCOND   DOUBLE PRECISION */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       DOUBLE PRECISION(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK SSIFA */
/*     BLAS SAXPY,SDOT,SSCAL,SASUM */
/*     FORTRAN DABS,DMAX1,IABS,DSIGN */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A USING ONLY UPPER HALF */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = sasum_(&j, &a[j * a_dim1 + 1], &c__1);
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] += (d__1 = a[i__ + j * a_dim1], abs(d__1));
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = z__[j];
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    ssifa_(&a[a_offset], lda, n, &kpvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE U*D*W = E */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L50: */
    }
    k = *n;
L60:
    if (k == 0) {
	goto L120;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L70;
    }
    t = z__[kps];
    z__[kps] = z__[kp];
    z__[kp] = t;
L70:
    if (z__[k] != 0.) {
	ek = d_sign(&ek, &z__[k]);
    }
    z__[k] += ek;
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 1) {
	goto L80;
    }
    if (z__[k - 1] != 0.) {
	ek = d_sign(&ek, &z__[k - 1]);
    }
    z__[k - 1] += ek;
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
	    c__1);
L80:
    if (ks == 2) {
	goto L100;
    }
    if ((d__1 = z__[k], abs(d__1)) <= (d__2 = a[k + k * a_dim1], abs(d__2))) {
	goto L90;
    }
    s = (d__1 = a[k + k * a_dim1], abs(d__1)) / (d__2 = z__[k], abs(d__2));
    sscal_(n, &s, &z__[1], &c__1);
    ek = s * ek;
L90:
    if (a[k + k * a_dim1] != 0.) {
	z__[k] /= a[k + k * a_dim1];
    }
    if (a[k + k * a_dim1] == 0.) {
	z__[k] = 1.;
    }
    goto L110;
L100:
    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    bk = z__[k] / a[k - 1 + k * a_dim1];
    bkm1 = z__[k - 1] / a[k - 1 + k * a_dim1];
    denom = ak * akm1 - 1.;
    z__[k] = (akm1 * bk - bkm1) / denom;
    z__[k - 1] = (ak * bkm1 - bk) / denom;
L110:
    k -= ks;
    goto L60;
L120:
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*     SOLVE TRANS(U)*Y = W */

    k = 1;
L130:
    if (k > *n) {
	goto L160;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L150;
    }
    i__1 = k - 1;
    z__[k] += sdot_(&i__1, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - 1;
	z__[k + 1] += sdot_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L140;
    }
    t = z__[k];
    z__[k] = z__[kp];
    z__[kp] = t;
L140:
L150:
    k += ks;
    goto L130;
L160:
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE U*D*V = Y */

    k = *n;
L170:
    if (k == 0) {
	goto L230;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == ks) {
	goto L190;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L180;
    }
    t = z__[kps];
    z__[kps] = z__[kp];
    z__[kp] = t;
L180:
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - ks;
	saxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
L190:
    if (ks == 2) {
	goto L210;
    }
    if ((d__1 = z__[k], abs(d__1)) <= (d__2 = a[k + k * a_dim1], abs(d__2))) {
	goto L200;
    }
    s = (d__1 = a[k + k * a_dim1], abs(d__1)) / (d__2 = z__[k], abs(d__2));
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;
L200:
    if (a[k + k * a_dim1] != 0.) {
	z__[k] /= a[k + k * a_dim1];
    }
    if (a[k + k * a_dim1] == 0.) {
	z__[k] = 1.;
    }
    goto L220;
L210:
    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    bk = z__[k] / a[k - 1 + k * a_dim1];
    bkm1 = z__[k - 1] / a[k - 1 + k * a_dim1];
    denom = ak * akm1 - 1.;
    z__[k] = (akm1 * bk - bkm1) / denom;
    z__[k - 1] = (ak * bkm1 - bk) / denom;
L220:
    k -= ks;
    goto L170;
L230:
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE TRANS(U)*Z = V */

    k = 1;
L240:
    if (k > *n) {
	goto L270;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L260;
    }
    i__1 = k - 1;
    z__[k] += sdot_(&i__1, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - 1;
	z__[k + 1] += sdot_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L250;
    }
    t = z__[k];
    z__[k] = z__[kp];
    z__[kp] = t;
L250:
L260:
    k += ks;
    goto L240;
L270:
/*     MAKE ZNORM = 1.0 */
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* ssico_ */

/* Subroutine */ int ssifa_(doublereal *a, integer *lda, integer *n, integer *
	kpvt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal t, ak, bk;
    static integer jj, km1, km2;
    static doublereal akm1, bkm1;
    static integer imax, jmax;
    static doublereal mulk;
    static logical swap;
    static doublereal alpha, denom;
    static integer kstep;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer imaxp1;
    static doublereal mulkm1, absakk;
    extern integer isamax_(integer *, doublereal *, integer *);
    static doublereal colmax, rowmax;


/*     SSIFA FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX BY ELIMINATION */
/*     WITH SYMMETRIC PIVOTING. */

/*     TO SOLVE  A*X = B , FOLLOW SSIFA BY DSISL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW SSIFA BY DSISL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW SSIFA BY DSIDI. */
/*     TO COMPUTE  INERTIA(A) , FOLLOW SSIFA BY DSIDI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW SSIFA BY DSIDI. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA,N) */
/*                THE SYMMETRIC MATRIX TO BE FACTORED. */
/*                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE */
/*                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS */
/*                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE, */
/*                     BUT IT DOES INDICATE THAT DSISL OR DSIDI MAY */
/*                     DIVIDE BY ZERO IF CALLED. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SSWAP,ISAMAX */
/*     FORTRAN DABS,DMAX1,DSQRT */

/*     INTERNAL VARIABLES */



/*     INITIALIZE */

/*     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;

    /* Function Body */
    alpha = (sqrt(17.) + 1.) / 8.;

    *info = 0;

/*     MAIN LOOP ON K, WHICH GOES FROM N TO 1. */

    k = *n;
L10:

/*        LEAVE THE LOOP IF K=0 OR K=1. */

/*     ...EXIT */
    if (k == 0) {
	goto L200;
    }
    if (k > 1) {
	goto L20;
    }
    kpvt[1] = 1;
    if (a[a_dim1 + 1] == 0.) {
	*info = 1;
    }
/*     ......EXIT */
    goto L200;
L20:

/*        THIS SECTION OF CODE DETERMINES THE KIND OF */
/*        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED, */
/*        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND */
/*        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS */
/*        REQUIRED. */

    km1 = k - 1;
    absakk = (d__1 = a[k + k * a_dim1], abs(d__1));

/*        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*        COLUMN K. */

    i__1 = k - 1;
    imax = isamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
    colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
    if (absakk < alpha * colmax) {
	goto L30;
    }
    kstep = 1;
    swap = FALSE_;
    goto L90;
L30:

/*           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*           ROW IMAX. */

    rowmax = 0.;
    imaxp1 = imax + 1;
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = rowmax, d__3 = (d__1 = a[imax + j * a_dim1], abs(d__1));
	rowmax = max(d__2,d__3);
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = isamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
    d__2 = rowmax, d__3 = (d__1 = a[jmax + imax * a_dim1], abs(d__1));
    rowmax = max(d__2,d__3);
L50:
    if ((d__1 = a[imax + imax * a_dim1], abs(d__1)) < alpha * rowmax) {
	goto L60;
    }
    kstep = 1;
    swap = TRUE_;
    goto L80;
L60:
    if (absakk < alpha * colmax * (colmax / rowmax)) {
	goto L70;
    }
    kstep = 1;
    swap = FALSE_;
    goto L80;
L70:
    kstep = 2;
    swap = imax != km1;
L80:
L90:
    if (max(absakk,colmax) != 0.) {
	goto L100;
    }

/*           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP. */

    kpvt[k] = k;
    *info = k;
    goto L190;
L100:
    if (kstep == 2) {
	goto L140;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (! swap) {
	goto L120;
    }

/*              PERFORM AN INTERCHANGE. */

    sswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	t = a[j + k * a_dim1];
	a[j + k * a_dim1] = a[imax + j * a_dim1];
	a[imax + j * a_dim1] = t;
/* L110: */
    }
L120:

/*           PERFORM THE ELIMINATION. */

    i__1 = km1;
    for (jj = 1; jj <= i__1; ++jj) {
	j = k - jj;
	mulk = -a[j + k * a_dim1] / a[k + k * a_dim1];
	t = mulk;
	saxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	a[j + k * a_dim1] = mulk;
/* L130: */
    }

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = k;
    if (swap) {
	kpvt[k] = imax;
    }
    goto L190;
L140:

/*           2 X 2 PIVOT BLOCK. */

    if (! swap) {
	goto L160;
    }

/*              PERFORM AN INTERCHANGE. */

    sswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[(k - 1) * a_dim1 + 1], &
	    c__1);
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	t = a[j + (k - 1) * a_dim1];
	a[j + (k - 1) * a_dim1] = a[imax + j * a_dim1];
	a[imax + j * a_dim1] = t;
/* L150: */
    }
    t = a[k - 1 + k * a_dim1];
    a[k - 1 + k * a_dim1] = a[imax + k * a_dim1];
    a[imax + k * a_dim1] = t;
L160:

/*           PERFORM THE ELIMINATION. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    denom = 1. - ak * akm1;
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	bk = a[j + k * a_dim1] / a[k - 1 + k * a_dim1];
	bkm1 = a[j + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
	mulk = (akm1 * bk - bkm1) / denom;
	mulkm1 = (ak * bkm1 - bk) / denom;
	t = mulk;
	saxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	t = mulkm1;
	saxpy_(&j, &t, &a[(k - 1) * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		c__1);
	a[j + k * a_dim1] = mulk;
	a[j + (k - 1) * a_dim1] = mulkm1;
/* L170: */
    }
L180:

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = 1 - k;
    if (swap) {
	kpvt[k] = -imax;
    }
    kpvt[k - 1] = kpvt[k];
L190:
    k -= kstep;
    goto L10;
L200:
    return 0;
} /* ssifa_ */

/* Subroutine */ int ssisl_(doublereal *a, integer *lda, integer *n, integer *
	kpvt, doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    static integer k;
    static doublereal ak, bk;
    static integer kp;
    static doublereal akm1, bkm1, temp;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal denom;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DSISL SOLVES THE DOUBLE PRECISION SYMMETRIC SYSTEM */
/*     A * X = B */
/*     USING THE FACTORS COMPUTED BY SSIFA. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA,N) */
/*                THE OUTPUT FROM SSIFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM SSIFA. */

/*        B       DOUBLE PRECISION(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO MAY OCCUR IF  DSICO  HAS SET RCOND .EQ. 0.0 */
/*        OR  SSIFA  HAS SET INFO .NE. 0  . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL SSIFA(A,LDA,N,KPVT,INFO) */
/*           IF (INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DSISL(A,LDA,N,KPVT,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SDOT */
/*     FORTRAN IABS */

/*     INTERNAL VARIABLES. */


/*     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND */
/*     D INVERSE TO B. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --b;

    /* Function Body */
    k = *n;
L10:
    if (k == 0) {
	goto L80;
    }
    if (kpvt[k] < 0) {
	goto L40;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L30;
    }
    kp = kpvt[k];
    if (kp == k) {
	goto L20;
    }

/*                 INTERCHANGE. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L20:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    saxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
L30:

/*           APPLY D INVERSE. */

    b[k] /= a[k + k * a_dim1];
    --k;
    goto L70;
L40:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 2) {
	goto L60;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k - 1) {
	goto L50;
    }

/*                 INTERCHANGE. */

    temp = b[k - 1];
    b[k - 1] = b[kp];
    b[kp] = temp;
L50:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 2;
    saxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    saxpy_(&i__1, &b[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
L60:

/*           APPLY D INVERSE. */

    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    bk = b[k] / a[k - 1 + k * a_dim1];
    bkm1 = b[k - 1] / a[k - 1 + k * a_dim1];
    denom = ak * akm1 - 1.;
    b[k] = (akm1 * bk - bkm1) / denom;
    b[k - 1] = (ak * bkm1 - bk) / denom;
    k += -2;
L70:
    goto L10;
L80:

/*     LOOP FORWARD APPLYING THE TRANSFORMATIONS. */

    k = 1;
L90:
    if (k > *n) {
	goto L160;
    }
    if (kpvt[k] < 0) {
	goto L120;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L110;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    b[k] += sdot_(&i__1, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    kp = kpvt[k];
    if (kp == k) {
	goto L100;
    }

/*                 INTERCHANGE. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L100:
L110:
    ++k;
    goto L150;
L120:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 1) {
	goto L140;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    b[k] += sdot_(&i__1, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    i__1 = k - 1;
    b[k + 1] += sdot_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L130;
    }

/*                 INTERCHANGE. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L130:
L140:
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* ssisl_ */

/* Subroutine */ int ssidi_(doublereal *a, integer *lda, integer *n, integer *
	kpvt, doublereal *det, integer *inert, doublereal *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal d__;
    static integer j, k;
    static doublereal t, ak;
    static integer jb, ks, km1;
    static doublereal ten, akp1, temp;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal akkp1;
    static logical nodet;
    static integer kstep;
    static logical noert, noinv;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);


/*     DSIDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE */
/*     OF A DOUBLE PRECISION SYMMETRIC MATRIX USING THE FACTORS FROM */
/*     SSIFA. */

/*     ON ENTRY */

/*        A       DOUBLE PRECISION(LDA,N) */
/*                THE OUTPUT FROM SSIFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY A. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX A. */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM SSIFA. */

/*        WORK    DOUBLE PRECISION(N) */
/*                WORK VECTOR.  CONTENTS DESTROYED. */

/*        JOB     INTEGER */
/*                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE */
/*                   IF  C .NE. 0, THE INVERSE IS COMPUTED, */
/*                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED, */
/*                   IF  A .NE. 0, THE INERTIA IS COMPUTED. */

/*                FOR EXAMPLE, JOB = 111  GIVES ALL THREE. */

/*     ON RETURN */

/*        VARIABLES NOT REQUESTED BY JOB ARE NOT USED. */

/*        A      CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF */
/*               THE ORIGINAL MATRIX.  THE STRICT LOWER TRIANGLE */
/*               IS NEVER REFERENCED. */

/*        DET    DOUBLE PRECISION(2) */
/*               DETERMINANT OF ORIGINAL MATRIX. */
/*               DETERMINANT = DET(1) * 10.0**DET(2) */
/*               WITH 1.0 .LE. DABS(DET(1)) .LT. 10.0 */
/*               OR DET(1) = 0.0. */

/*        INERT  INTEGER(3) */
/*               THE INERTIA OF THE ORIGINAL MATRIX. */
/*               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES. */
/*               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES. */
/*               INERT(3)  =  NUMBER OF ZERO EIGENVALUES. */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO MAY OCCUR IF THE INVERSE IS REQUESTED */
/*        AND  DSICO  HAS SET RCOND .EQ. 0.0 */
/*        OR  SSIFA  HAS SET  INFO .NE. 0 . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SCOPY,SDOT,SSWAP */
/*     FORTRAN DABS,IABS,MOD */

/*     INTERNAL VARIABLES. */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --det;
    --inert;
    --work;

    /* Function Body */
    noinv = *job % 10 == 0;
    nodet = *job % 100 / 10 == 0;
    noert = *job % 1000 / 100 == 0;

    if (nodet && noert) {
	goto L140;
    }
    if (noert) {
	goto L10;
    }
    inert[1] = 0;
    inert[2] = 0;
    inert[3] = 0;
L10:
    if (nodet) {
	goto L20;
    }
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
L20:
    t = 0.;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	d__ = a[k + k * a_dim1];

/*           CHECK IF 1 BY 1 */

	if (kpvt[k] > 0) {
	    goto L50;
	}

/*              2 BY 2 BLOCK */
/*              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S) */
/*                      (S  C) */
/*              TO AVOID UNDERFLOW/OVERFLOW TROUBLES. */
/*              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG. */

	if (t != 0.) {
	    goto L30;
	}
	t = (d__1 = a[k + (k + 1) * a_dim1], abs(d__1));
	d__ = d__ / t * a[k + 1 + (k + 1) * a_dim1] - t;
	goto L40;
L30:
	d__ = t;
	t = 0.;
L40:
L50:

	if (noert) {
	    goto L60;
	}
	if (d__ > 0.) {
	    ++inert[1];
	}
	if (d__ < 0.) {
	    ++inert[2];
	}
	if (d__ == 0.) {
	    ++inert[3];
	}
L60:

	if (nodet) {
	    goto L120;
	}
	det[1] = d__ * det[1];
	if (det[1] == 0.) {
	    goto L110;
	}
L70:
	if (abs(det[1]) >= 1.) {
	    goto L80;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L70;
L80:
L90:
	if (abs(det[1]) < ten) {
	    goto L100;
	}
	det[1] /= ten;
	det[2] += 1.;
	goto L90;
L100:
L110:
L120:
/* L130: */
	;
    }
L140:

/*     COMPUTE INVERSE(A) */

    if (noinv) {
	goto L270;
    }
    k = 1;
L150:
    if (k > *n) {
	goto L260;
    }
    km1 = k - 1;
    if (kpvt[k] < 0) {
	goto L180;
    }

/*              1 BY 1 */

    a[k + k * a_dim1] = 1. / a[k + k * a_dim1];
    if (km1 < 1) {
	goto L170;
    }
    scopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + k * a_dim1] = sdot_(&j, &a[j * a_dim1 + 1], &c__1, &work[1], &
		c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L160: */
    }
    a[k + k * a_dim1] += sdot_(&km1, &work[1], &c__1, &a[k * a_dim1 + 1], &
	    c__1);
L170:
    kstep = 1;
    goto L220;
L180:

/*              2 BY 2 */

    t = (d__1 = a[k + (k + 1) * a_dim1], abs(d__1));
    ak = a[k + k * a_dim1] / t;
    akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
    akkp1 = a[k + (k + 1) * a_dim1] / t;
    d__ = t * (ak * akp1 - 1.);
    a[k + k * a_dim1] = akp1 / d__;
    a[k + 1 + (k + 1) * a_dim1] = ak / d__;
    a[k + (k + 1) * a_dim1] = -akkp1 / d__;
    if (km1 < 1) {
	goto L210;
    }
    scopy_(&km1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + (k + 1) * a_dim1] = sdot_(&j, &a[j * a_dim1 + 1], &c__1, &work[
		1], &c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[(k + 1) * 
		a_dim1 + 1], &c__1);
/* L190: */
    }
    a[k + 1 + (k + 1) * a_dim1] += sdot_(&km1, &work[1], &c__1, &a[(k + 1) * 
	    a_dim1 + 1], &c__1);
    a[k + (k + 1) * a_dim1] += sdot_(&km1, &a[k * a_dim1 + 1], &c__1, &a[(k + 
	    1) * a_dim1 + 1], &c__1);
    scopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + k * a_dim1] = sdot_(&j, &a[j * a_dim1 + 1], &c__1, &work[1], &
		c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L200: */
    }
    a[k + k * a_dim1] += sdot_(&km1, &work[1], &c__1, &a[k * a_dim1 + 1], &
	    c__1);
L210:
    kstep = 2;
L220:

/*           SWAP */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L250;
    }
    sswap_(&ks, &a[ks * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	temp = a[j + k * a_dim1];
	a[j + k * a_dim1] = a[ks + j * a_dim1];
	a[ks + j * a_dim1] = temp;
/* L230: */
    }
    if (kstep == 1) {
	goto L240;
    }
    temp = a[ks + (k + 1) * a_dim1];
    a[ks + (k + 1) * a_dim1] = a[k + (k + 1) * a_dim1];
    a[k + (k + 1) * a_dim1] = temp;
L240:
L250:
    k += kstep;
    goto L150;
L260:
L270:
    return 0;
} /* ssidi_ */

/* Subroutine */ int sspco_(doublereal *ap, integer *n, integer *kpvt, 
	doublereal *rcond, doublereal *z__)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s, t;
    static integer j1;
    static doublereal ak, bk, ek;
    static integer ij, ik, kk, kp, ks, jm1, kps;
    static doublereal akm1, bkm1;
    static integer ikm1, km1k, ikp1, info;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer km1km1;
    static doublereal denom;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal anorm;
    extern /* Subroutine */ int sspfa_(doublereal *, integer *, integer *, 
	    integer *);
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DSPCO FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX STORED IN */
/*     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING AND ESTIMATES */
/*     THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, SSPFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW DSPCO BY DSPSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW DSPCO BY DSPSL. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW DSPCO BY DSPDI. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW DSPCO BY DSPDI. */
/*     TO COMPUTE  INERTIA(A), FOLLOW DSPCO BY DSPDI. */

/*     ON ENTRY */

/*        AP      DOUBLE PRECISION (N*(N+1)/2) */
/*                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE */
/*                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY */
/*                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 . */
/*                SEE COMMENTS BELOW FOR DETAILS. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     OUTPUT */

/*        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT STORED IN PACKED FORM. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE */
/*                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        RCOND   DOUBLE PRECISION */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       DOUBLE PRECISION(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     PACKED STORAGE */

/*          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER */
/*          TRIANGLE OF A SYMMETRIC MATRIX. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K) = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK SSPFA */
/*     BLAS SAXPY,SDOT,SSCAL,SASUM */
/*     FORTRAN DABS,DMAX1,IABS,DSIGN */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A USING ONLY UPPER HALF */

    /* Parameter adjustments */
    --z__;
    --kpvt;
    --ap;

    /* Function Body */
    j1 = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = sasum_(&j, &ap[j1], &c__1);
	ij = j1;
	j1 += j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] += (d__1 = ap[ij], abs(d__1));
	    ++ij;
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = z__[j];
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    sspfa_(&ap[1], n, &kpvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE U*D*W = E */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L50: */
    }
    k = *n;
    ik = *n * (*n - 1) / 2;
L60:
    if (k == 0) {
	goto L120;
    }
    kk = ik + k;
    ikm1 = ik - (k - 1);
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L70;
    }
    t = z__[kps];
    z__[kps] = z__[kp];
    z__[kp] = t;
L70:
    if (z__[k] != 0.) {
	ek = d_sign(&ek, &z__[k]);
    }
    z__[k] += ek;
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k], &ap[ik + 1], &c__1, &z__[1], &c__1);
    if (ks == 1) {
	goto L80;
    }
    if (z__[k - 1] != 0.) {
	ek = d_sign(&ek, &z__[k - 1]);
    }
    z__[k - 1] += ek;
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k - 1], &ap[ikm1 + 1], &c__1, &z__[1], &c__1);
L80:
    if (ks == 2) {
	goto L100;
    }
    if ((d__1 = z__[k], abs(d__1)) <= (d__2 = ap[kk], abs(d__2))) {
	goto L90;
    }
    s = (d__1 = ap[kk], abs(d__1)) / (d__2 = z__[k], abs(d__2));
    sscal_(n, &s, &z__[1], &c__1);
    ek = s * ek;
L90:
    if (ap[kk] != 0.) {
	z__[k] /= ap[kk];
    }
    if (ap[kk] == 0.) {
	z__[k] = 1.;
    }
    goto L110;
L100:
    km1k = ik + k - 1;
    km1km1 = ikm1 + k - 1;
    ak = ap[kk] / ap[km1k];
    akm1 = ap[km1km1] / ap[km1k];
    bk = z__[k] / ap[km1k];
    bkm1 = z__[k - 1] / ap[km1k];
    denom = ak * akm1 - 1.;
    z__[k] = (akm1 * bk - bkm1) / denom;
    z__[k - 1] = (ak * bkm1 - bk) / denom;
L110:
    k -= ks;
    ik -= k;
    if (ks == 2) {
	ik -= k + 1;
    }
    goto L60;
L120:
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*     SOLVE TRANS(U)*Y = W */

    k = 1;
    ik = 0;
L130:
    if (k > *n) {
	goto L160;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L150;
    }
    i__1 = k - 1;
    z__[k] += sdot_(&i__1, &ap[ik + 1], &c__1, &z__[1], &c__1);
    ikp1 = ik + k;
    if (ks == 2) {
	i__1 = k - 1;
	z__[k + 1] += sdot_(&i__1, &ap[ikp1 + 1], &c__1, &z__[1], &c__1);
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L140;
    }
    t = z__[k];
    z__[k] = z__[kp];
    z__[kp] = t;
L140:
L150:
    ik += k;
    if (ks == 2) {
	ik += k + 1;
    }
    k += ks;
    goto L130;
L160:
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE U*D*V = Y */

    k = *n;
    ik = *n * (*n - 1) / 2;
L170:
    if (k == 0) {
	goto L230;
    }
    kk = ik + k;
    ikm1 = ik - (k - 1);
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == ks) {
	goto L190;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L180;
    }
    t = z__[kps];
    z__[kps] = z__[kp];
    z__[kp] = t;
L180:
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k], &ap[ik + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - ks;
	saxpy_(&i__1, &z__[k - 1], &ap[ikm1 + 1], &c__1, &z__[1], &c__1);
    }
L190:
    if (ks == 2) {
	goto L210;
    }
    if ((d__1 = z__[k], abs(d__1)) <= (d__2 = ap[kk], abs(d__2))) {
	goto L200;
    }
    s = (d__1 = ap[kk], abs(d__1)) / (d__2 = z__[k], abs(d__2));
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;
L200:
    if (ap[kk] != 0.) {
	z__[k] /= ap[kk];
    }
    if (ap[kk] == 0.) {
	z__[k] = 1.;
    }
    goto L220;
L210:
    km1k = ik + k - 1;
    km1km1 = ikm1 + k - 1;
    ak = ap[kk] / ap[km1k];
    akm1 = ap[km1km1] / ap[km1k];
    bk = z__[k] / ap[km1k];
    bkm1 = z__[k - 1] / ap[km1k];
    denom = ak * akm1 - 1.;
    z__[k] = (akm1 * bk - bkm1) / denom;
    z__[k - 1] = (ak * bkm1 - bk) / denom;
L220:
    k -= ks;
    ik -= k;
    if (ks == 2) {
	ik -= k + 1;
    }
    goto L170;
L230:
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE TRANS(U)*Z = V */

    k = 1;
    ik = 0;
L240:
    if (k > *n) {
	goto L270;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L260;
    }
    i__1 = k - 1;
    z__[k] += sdot_(&i__1, &ap[ik + 1], &c__1, &z__[1], &c__1);
    ikp1 = ik + k;
    if (ks == 2) {
	i__1 = k - 1;
	z__[k + 1] += sdot_(&i__1, &ap[ikp1 + 1], &c__1, &z__[1], &c__1);
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L250;
    }
    t = z__[k];
    z__[k] = z__[kp];
    z__[kp] = t;
L250:
L260:
    ik += k;
    if (ks == 2) {
	ik += k + 1;
    }
    k += ks;
    goto L240;
L270:
/*     MAKE ZNORM = 1.0 */
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* sspco_ */

/* Subroutine */ int sspfa_(doublereal *ap, integer *n, integer *kpvt, 
	integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal t, ak, bk;
    static integer ij, ik, jj, im, jk, kk, km1, km2, ijj, imj, imk;
    static doublereal akm1, bkm1;
    static integer ikm1, jkm1, km1k, imim, jmim, imax, jmax;
    static doublereal mulk;
    static logical swap;
    static doublereal alpha;
    static integer km1km1;
    static doublereal denom;
    static integer kstep;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer imaxp1;
    static doublereal mulkm1, absakk;
    extern integer isamax_(integer *, doublereal *, integer *);
    static doublereal colmax, rowmax;


/*     SSPFA FACTORS A DOUBLE PRECISION SYMMETRIC MATRIX STORED IN */
/*     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING. */

/*     TO SOLVE  A*X = B , FOLLOW SSPFA BY DSPSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW SSPFA BY DSPSL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW SSPFA BY DSPDI. */
/*     TO COMPUTE  INERTIA(A) , FOLLOW SSPFA BY DSPDI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW SSPFA BY DSPDI. */

/*     ON ENTRY */

/*        AP      DOUBLE PRECISION (N*(N+1)/2) */
/*                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE */
/*                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY */
/*                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 . */
/*                SEE COMMENTS BELOW FOR DETAILS. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     OUTPUT */

/*        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT STORED IN PACKED FORM. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE */
/*                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS */
/*                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE, */
/*                     BUT IT DOES INDICATE THAT DSPSL OR DSPDI MAY */
/*                     DIVIDE BY ZERO IF CALLED. */

/*     PACKED STORAGE */

/*          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER */
/*          TRIANGLE OF A SYMMETRIC MATRIX. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K)  = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SSWAP,ISAMAX */
/*     FORTRAN DABS,DMAX1,DSQRT */

/*     INTERNAL VARIABLES */



/*     INITIALIZE */

/*     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE. */
    /* Parameter adjustments */
    --kpvt;
    --ap;

    /* Function Body */
    alpha = (sqrt(17.) + 1.) / 8.;

    *info = 0;

/*     MAIN LOOP ON K, WHICH GOES FROM N TO 1. */

    k = *n;
    ik = *n * (*n - 1) / 2;
L10:

/*        LEAVE THE LOOP IF K=0 OR K=1. */

/*     ...EXIT */
    if (k == 0) {
	goto L200;
    }
    if (k > 1) {
	goto L20;
    }
    kpvt[1] = 1;
    if (ap[1] == 0.) {
	*info = 1;
    }
/*     ......EXIT */
    goto L200;
L20:

/*        THIS SECTION OF CODE DETERMINES THE KIND OF */
/*        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED, */
/*        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND */
/*        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS */
/*        REQUIRED. */

    km1 = k - 1;
    kk = ik + k;
    absakk = (d__1 = ap[kk], abs(d__1));

/*        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*        COLUMN K. */

    i__1 = k - 1;
    imax = isamax_(&i__1, &ap[ik + 1], &c__1);
    imk = ik + imax;
    colmax = (d__1 = ap[imk], abs(d__1));
    if (absakk < alpha * colmax) {
	goto L30;
    }
    kstep = 1;
    swap = FALSE_;
    goto L90;
L30:

/*           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*           ROW IMAX. */

    rowmax = 0.;
    imaxp1 = imax + 1;
    im = imax * (imax - 1) / 2;
    imj = im + (imax << 1);
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = rowmax, d__3 = (d__1 = ap[imj], abs(d__1));
	rowmax = max(d__2,d__3);
	imj += j;
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = isamax_(&i__1, &ap[im + 1], &c__1);
    jmim = jmax + im;
/* Computing MAX */
    d__2 = rowmax, d__3 = (d__1 = ap[jmim], abs(d__1));
    rowmax = max(d__2,d__3);
L50:
    imim = imax + im;
    if ((d__1 = ap[imim], abs(d__1)) < alpha * rowmax) {
	goto L60;
    }
    kstep = 1;
    swap = TRUE_;
    goto L80;
L60:
    if (absakk < alpha * colmax * (colmax / rowmax)) {
	goto L70;
    }
    kstep = 1;
    swap = FALSE_;
    goto L80;
L70:
    kstep = 2;
    swap = imax != km1;
L80:
L90:
    if (max(absakk,colmax) != 0.) {
	goto L100;
    }

/*           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP. */

    kpvt[k] = k;
    *info = k;
    goto L190;
L100:
    if (kstep == 2) {
	goto L140;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (! swap) {
	goto L120;
    }

/*              PERFORM AN INTERCHANGE. */

    sswap_(&imax, &ap[im + 1], &c__1, &ap[ik + 1], &c__1);
    imj = ik + imax;
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	jk = ik + j;
	t = ap[jk];
	ap[jk] = ap[imj];
	ap[imj] = t;
	imj -= j - 1;
/* L110: */
    }
L120:

/*           PERFORM THE ELIMINATION. */

    ij = ik - (k - 1);
    i__1 = km1;
    for (jj = 1; jj <= i__1; ++jj) {
	j = k - jj;
	jk = ik + j;
	mulk = -ap[jk] / ap[kk];
	t = mulk;
	saxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	ijj = ij + j;
	ap[jk] = mulk;
	ij -= j - 1;
/* L130: */
    }

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = k;
    if (swap) {
	kpvt[k] = imax;
    }
    goto L190;
L140:

/*           2 X 2 PIVOT BLOCK. */

    km1k = ik + k - 1;
    ikm1 = ik - (k - 1);
    if (! swap) {
	goto L160;
    }

/*              PERFORM AN INTERCHANGE. */

    sswap_(&imax, &ap[im + 1], &c__1, &ap[ikm1 + 1], &c__1);
    imj = ikm1 + imax;
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	jkm1 = ikm1 + j;
	t = ap[jkm1];
	ap[jkm1] = ap[imj];
	ap[imj] = t;
	imj -= j - 1;
/* L150: */
    }
    t = ap[km1k];
    ap[km1k] = ap[imk];
    ap[imk] = t;
L160:

/*           PERFORM THE ELIMINATION. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    ak = ap[kk] / ap[km1k];
    km1km1 = ikm1 + k - 1;
    akm1 = ap[km1km1] / ap[km1k];
    denom = 1. - ak * akm1;
    ij = ik - (k - 1) - (k - 2);
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	jk = ik + j;
	bk = ap[jk] / ap[km1k];
	jkm1 = ikm1 + j;
	bkm1 = ap[jkm1] / ap[km1k];
	mulk = (akm1 * bk - bkm1) / denom;
	mulkm1 = (ak * bkm1 - bk) / denom;
	t = mulk;
	saxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	t = mulkm1;
	saxpy_(&j, &t, &ap[ikm1 + 1], &c__1, &ap[ij + 1], &c__1);
	ap[jk] = mulk;
	ap[jkm1] = mulkm1;
	ijj = ij + j;
	ij -= j - 1;
/* L170: */
    }
L180:

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = 1 - k;
    if (swap) {
	kpvt[k] = -imax;
    }
    kpvt[k - 1] = kpvt[k];
L190:
    ik -= k - 1;
    if (kstep == 2) {
	ik -= k - 2;
    }
    k -= kstep;
    goto L10;
L200:
    return 0;
} /* sspfa_ */

/* Subroutine */ int sspsl_(doublereal *ap, integer *n, integer *kpvt, 
	doublereal *b)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    static doublereal ak, bk;
    static integer ik, kk, kp;
    static doublereal akm1, bkm1;
    static integer ikm1, km1k, ikp1;
    static doublereal temp;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer km1km1;
    static doublereal denom;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DSISL SOLVES THE DOUBLE PRECISION SYMMETRIC SYSTEM */
/*     A * X = B */
/*     USING THE FACTORS COMPUTED BY SSPFA. */

/*     ON ENTRY */

/*        AP      DOUBLE PRECISION(N*(N+1)/2) */
/*                THE OUTPUT FROM SSPFA. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM SSPFA. */

/*        B       DOUBLE PRECISION(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO MAY OCCUR IF  DSPCO  HAS SET RCOND .EQ. 0.0 */
/*        OR  SSPFA  HAS SET INFO .NE. 0  . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL SSPFA(AP,N,KPVT,INFO) */
/*           IF (INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL DSPSL(AP,N,KPVT,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SDOT */
/*     FORTRAN IABS */

/*     INTERNAL VARIABLES. */


/*     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND */
/*     D INVERSE TO B. */

    /* Parameter adjustments */
    --b;
    --kpvt;
    --ap;

    /* Function Body */
    k = *n;
    ik = *n * (*n - 1) / 2;
L10:
    if (k == 0) {
	goto L80;
    }
    kk = ik + k;
    if (kpvt[k] < 0) {
	goto L40;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L30;
    }
    kp = kpvt[k];
    if (kp == k) {
	goto L20;
    }

/*                 INTERCHANGE. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L20:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    saxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
L30:

/*           APPLY D INVERSE. */

    b[k] /= ap[kk];
    --k;
    ik -= k;
    goto L70;
L40:

/*           2 X 2 PIVOT BLOCK. */

    ikm1 = ik - (k - 1);
    if (k == 2) {
	goto L60;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k - 1) {
	goto L50;
    }

/*                 INTERCHANGE. */

    temp = b[k - 1];
    b[k - 1] = b[kp];
    b[kp] = temp;
L50:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 2;
    saxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    saxpy_(&i__1, &b[k - 1], &ap[ikm1 + 1], &c__1, &b[1], &c__1);
L60:

/*           APPLY D INVERSE. */

    km1k = ik + k - 1;
    kk = ik + k;
    ak = ap[kk] / ap[km1k];
    km1km1 = ikm1 + k - 1;
    akm1 = ap[km1km1] / ap[km1k];
    bk = b[k] / ap[km1k];
    bkm1 = b[k - 1] / ap[km1k];
    denom = ak * akm1 - 1.;
    b[k] = (akm1 * bk - bkm1) / denom;
    b[k - 1] = (ak * bkm1 - bk) / denom;
    k += -2;
    ik = ik - (k + 1) - k;
L70:
    goto L10;
L80:

/*     LOOP FORWARD APPLYING THE TRANSFORMATIONS. */

    k = 1;
    ik = 0;
L90:
    if (k > *n) {
	goto L160;
    }
    if (kpvt[k] < 0) {
	goto L120;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L110;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    b[k] += sdot_(&i__1, &ap[ik + 1], &c__1, &b[1], &c__1);
    kp = kpvt[k];
    if (kp == k) {
	goto L100;
    }

/*                 INTERCHANGE. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L100:
L110:
    ik += k;
    ++k;
    goto L150;
L120:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 1) {
	goto L140;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    b[k] += sdot_(&i__1, &ap[ik + 1], &c__1, &b[1], &c__1);
    ikp1 = ik + k;
    i__1 = k - 1;
    b[k + 1] += sdot_(&i__1, &ap[ikp1 + 1], &c__1, &b[1], &c__1);
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L130;
    }

/*                 INTERCHANGE. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L130:
L140:
    ik = ik + k + k + 1;
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* sspsl_ */

/* Subroutine */ int sspdi_(doublereal *ap, integer *n, integer *kpvt, 
	doublereal *det, integer *inert, doublereal *work, integer *job)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static doublereal d__;
    static integer j, k;
    static doublereal t, ak;
    static integer jb, ij, ik, jk, kk, ks, km1;
    static doublereal ten;
    static integer iks, ksj;
    static doublereal akp1;
    static integer ikp1, jkp1, kkp1;
    static doublereal temp;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal akkp1;
    static integer kskp1;
    static logical nodet;
    static integer kstep;
    static logical noert, noinv;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);


/*     DSPDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE */
/*     OF A DOUBLE PRECISION SYMMETRIC MATRIX USING THE FACTORS FROM */
/*     SSPFA, WHERE THE MATRIX IS STORED IN PACKED FORM. */

/*     ON ENTRY */

/*        AP      DOUBLE PRECISION (N*(N+1)/2) */
/*                THE OUTPUT FROM SSPFA. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX A. */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM SSPFA. */

/*        WORK    DOUBLE PRECISION(N) */
/*                WORK VECTOR.  CONTENTS IGNORED. */

/*        JOB     INTEGER */
/*                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE */
/*                   IF  C .NE. 0, THE INVERSE IS COMPUTED, */
/*                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED, */
/*                   IF  A .NE. 0, THE INERTIA IS COMPUTED. */

/*                FOR EXAMPLE, JOB = 111  GIVES ALL THREE. */

/*     ON RETURN */

/*        VARIABLES NOT REQUESTED BY JOB ARE NOT USED. */

/*        AP     CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF */
/*               THE ORIGINAL MATRIX, STORED IN PACKED FORM. */
/*               THE COLUMNS OF THE UPPER TRIANGLE ARE STORED */
/*               SEQUENTIALLY IN A ONE-DIMENSIONAL ARRAY. */

/*        DET    DOUBLE PRECISION(2) */
/*               DETERMINANT OF ORIGINAL MATRIX. */
/*               DETERMINANT = DET(1) * 10.0**DET(2) */
/*               WITH 1.0 .LE. DABS(DET(1)) .LT. 10.0 */
/*               OR DET(1) = 0.0. */

/*        INERT  INTEGER(3) */
/*               THE INERTIA OF THE ORIGINAL MATRIX. */
/*               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES. */
/*               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES. */
/*               INERT(3)  =  NUMBER OF ZERO EIGENVALUES. */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INVERSE IS REQUESTED */
/*        AND  DSPCO  HAS SET RCOND .EQ. 0.0 */
/*        OR  SSPFA  HAS SET  INFO .NE. 0 . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SCOPY,SDOT,SSWAP */
/*     FORTRAN DABS,IABS,MOD */

/*     INTERNAL VARIABLES. */


    /* Parameter adjustments */
    --work;
    --inert;
    --det;
    --kpvt;
    --ap;

    /* Function Body */
    noinv = *job % 10 == 0;
    nodet = *job % 100 / 10 == 0;
    noert = *job % 1000 / 100 == 0;

    if (nodet && noert) {
	goto L140;
    }
    if (noert) {
	goto L10;
    }
    inert[1] = 0;
    inert[2] = 0;
    inert[3] = 0;
L10:
    if (nodet) {
	goto L20;
    }
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
L20:
    t = 0.;
    ik = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk = ik + k;
	d__ = ap[kk];

/*           CHECK IF 1 BY 1 */

	if (kpvt[k] > 0) {
	    goto L50;
	}

/*              2 BY 2 BLOCK */
/*              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = DABS(S) */
/*                      (S  C) */
/*              TO AVOID UNDERFLOW/OVERFLOW TROUBLES. */
/*              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG. */

	if (t != 0.) {
	    goto L30;
	}
	ikp1 = ik + k;
	kkp1 = ikp1 + k;
	t = (d__1 = ap[kkp1], abs(d__1));
	d__ = d__ / t * ap[kkp1 + 1] - t;
	goto L40;
L30:
	d__ = t;
	t = 0.;
L40:
L50:

	if (noert) {
	    goto L60;
	}
	if (d__ > 0.) {
	    ++inert[1];
	}
	if (d__ < 0.) {
	    ++inert[2];
	}
	if (d__ == 0.) {
	    ++inert[3];
	}
L60:

	if (nodet) {
	    goto L120;
	}
	det[1] = d__ * det[1];
	if (det[1] == 0.) {
	    goto L110;
	}
L70:
	if (abs(det[1]) >= 1.) {
	    goto L80;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L70;
L80:
L90:
	if (abs(det[1]) < ten) {
	    goto L100;
	}
	det[1] /= ten;
	det[2] += 1.;
	goto L90;
L100:
L110:
L120:
	ik += k;
/* L130: */
    }
L140:

/*     COMPUTE INVERSE(A) */

    if (noinv) {
	goto L270;
    }
    k = 1;
    ik = 0;
L150:
    if (k > *n) {
	goto L260;
    }
    km1 = k - 1;
    kk = ik + k;
    ikp1 = ik + k;
    kkp1 = ikp1 + k;
    if (kpvt[k] < 0) {
	goto L180;
    }

/*              1 BY 1 */

    ap[kk] = 1. / ap[kk];
    if (km1 < 1) {
	goto L170;
    }
    scopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	ap[jk] = sdot_(&j, &ap[ij + 1], &c__1, &work[1], &c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L160: */
    }
    ap[kk] += sdot_(&km1, &work[1], &c__1, &ap[ik + 1], &c__1);
L170:
    kstep = 1;
    goto L220;
L180:

/*              2 BY 2 */

    t = (d__1 = ap[kkp1], abs(d__1));
    ak = ap[kk] / t;
    akp1 = ap[kkp1 + 1] / t;
    akkp1 = ap[kkp1] / t;
    d__ = t * (ak * akp1 - 1.);
    ap[kk] = akp1 / d__;
    ap[kkp1 + 1] = ak / d__;
    ap[kkp1] = -akkp1 / d__;
    if (km1 < 1) {
	goto L210;
    }
    scopy_(&km1, &ap[ikp1 + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jkp1 = ikp1 + j;
	ap[jkp1] = sdot_(&j, &ap[ij + 1], &c__1, &work[1], &c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ikp1 + 1], &c__1);
	ij += j;
/* L190: */
    }
    ap[kkp1 + 1] += sdot_(&km1, &work[1], &c__1, &ap[ikp1 + 1], &c__1);
    ap[kkp1] += sdot_(&km1, &ap[ik + 1], &c__1, &ap[ikp1 + 1], &c__1);
    scopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	ap[jk] = sdot_(&j, &ap[ij + 1], &c__1, &work[1], &c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L200: */
    }
    ap[kk] += sdot_(&km1, &work[1], &c__1, &ap[ik + 1], &c__1);
L210:
    kstep = 2;
L220:

/*           SWAP */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L250;
    }
    iks = ks * (ks - 1) / 2;
    sswap_(&ks, &ap[iks + 1], &c__1, &ap[ik + 1], &c__1);
    ksj = ik + ks;
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	jk = ik + j;
	temp = ap[jk];
	ap[jk] = ap[ksj];
	ap[ksj] = temp;
	ksj -= j - 1;
/* L230: */
    }
    if (kstep == 1) {
	goto L240;
    }
    kskp1 = ikp1 + ks;
    temp = ap[kskp1];
    ap[kskp1] = ap[kkp1];
    ap[kkp1] = temp;
L240:
L250:
    ik += k;
    if (kstep == 2) {
	ik = ik + k + 1;
    }
    k += kstep;
    goto L150;
L260:
L270:
    return 0;
} /* sspdi_ */

/* Subroutine */ int strco_(doublereal *t, integer *ldt, integer *n, 
	doublereal *rcond, doublereal *z__, integer *job)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer j, k, l;
    static doublereal s, w;
    static integer i1, j1, j2;
    static doublereal ek;
    static integer kk;
    static doublereal sm, wk, wkm;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern doublereal sasum_(integer *, doublereal *, integer *);
    static logical lower;
    static doublereal tnorm, ynorm;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     DTRCO ESTIMATES THE CONDITION OF A DOUBLE PRECISION TRIANGULAR */
/*     MATRIX. */

/*     ON ENTRY */

/*        T       DOUBLE PRECISION(LDT,N) */
/*                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO */
/*                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND */
/*                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE */
/*                USED TO STORE OTHER INFORMATION. */

/*        LDT     INTEGER */
/*                LDT IS THE LEADING DIMENSION OF THE ARRAY T. */

/*        N       INTEGER */
/*                N IS THE ORDER OF THE SYSTEM. */

/*        JOB     INTEGER */
/*                = 0         T  IS LOWER TRIANGULAR. */
/*                = NONZERO   T  IS UPPER TRIANGULAR. */

/*     ON RETURN */

/*        RCOND   DOUBLE PRECISION */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  T . */
/*                FOR THE SYSTEM  T*X = B , RELATIVE PERTURBATIONS */
/*                IN  T  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  T  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       DOUBLE PRECISION(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  T  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SSCAL,SASUM */
/*     FORTRAN DABS,DMAX1,DSIGN */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --z__;

    /* Function Body */
    lower = *job == 0;

/*     COMPUTE 1-NORM OF T */

    tnorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = j;
	if (lower) {
	    l = *n + 1 - j;
	}
	i1 = 1;
	if (lower) {
	    i1 = j;
	}
/* Computing MAX */
	d__1 = tnorm, d__2 = sasum_(&l, &t[i1 + j * t_dim1], &c__1);
	tnorm = max(d__1,d__2);
/* L10: */
    }

/*     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E . */
/*     TRANS(T)  IS THE TRANSPOSE OF T . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF Y . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE TRANS(T)*Y = E */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L20: */
    }
    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = kk;
	if (lower) {
	    k = *n + 1 - kk;
	}
	if (z__[k] != 0.) {
	    d__1 = -z__[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z__[k], abs(d__1)) <= (d__2 = t[k + k * t_dim1], abs(
		d__2))) {
	    goto L30;
	}
	s = (d__1 = t[k + k * t_dim1], abs(d__1)) / (d__2 = ek - z__[k], abs(
		d__2));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = abs(wk);
	sm = abs(wkm);
	if (t[k + k * t_dim1] == 0.) {
	    goto L40;
	}
	wk /= t[k + k * t_dim1];
	wkm /= t[k + k * t_dim1];
	goto L50;
L40:
	wk = 1.;
	wkm = 1.;
L50:
	if (kk == *n) {
	    goto L90;
	}
	j1 = k + 1;
	if (lower) {
	    j1 = 1;
	}
	j2 = *n;
	if (lower) {
	    j2 = k - 1;
	}
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    sm += (d__1 = z__[j] + wkm * t[k + j * t_dim1], abs(d__1));
	    z__[j] += wk * t[k + j * t_dim1];
	    s += (d__1 = z__[j], abs(d__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	w = wkm - wk;
	wk = wkm;
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    z__[j] += w * t[k + j * t_dim1];
/* L70: */
	}
L80:
L90:
	z__[k] = wk;
/* L100: */
    }
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE T*Z = Y */

    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = *n + 1 - kk;
	if (lower) {
	    k = kk;
	}
	if ((d__1 = z__[k], abs(d__1)) <= (d__2 = t[k + k * t_dim1], abs(d__2)
		)) {
	    goto L110;
	}
	s = (d__1 = t[k + k * t_dim1], abs(d__1)) / (d__2 = z__[k], abs(d__2))
		;
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L110:
	if (t[k + k * t_dim1] != 0.) {
	    z__[k] /= t[k + k * t_dim1];
	}
	if (t[k + k * t_dim1] == 0.) {
	    z__[k] = 1.;
	}
	i1 = 1;
	if (lower) {
	    i1 = k + 1;
	}
	if (kk >= *n) {
	    goto L120;
	}
	w = -z__[k];
	i__2 = *n - kk;
	saxpy_(&i__2, &w, &t[i1 + k * t_dim1], &c__1, &z__[i1], &c__1);
L120:
/* L130: */
	;
    }
/*     MAKE ZNORM = 1.0 */
    s = 1. / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (tnorm != 0.) {
	*rcond = ynorm / tnorm;
    }
    if (tnorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* strco_ */

/* Subroutine */ int strsl_(doublereal *t, integer *ldt, integer *n, 
	doublereal *b, integer *job, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2;

    /* Local variables */
    static integer j, jj, case__;
    static doublereal temp;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);



/*     DTRSL SOLVES SYSTEMS OF THE FORM */

/*                   T * X = B */
/*     OR */
/*                   TRANS(T) * X = B */

/*     WHERE T IS A TRIANGULAR MATRIX OF ORDER N. HERE TRANS(T) */
/*     DENOTES THE TRANSPOSE OF THE MATRIX T. */

/*     ON ENTRY */

/*         T         DOUBLE PRECISION(LDT,N) */
/*                   T CONTAINS THE MATRIX OF THE SYSTEM. THE ZERO */
/*                   ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND */
/*                   THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE */
/*                   USED TO STORE OTHER INFORMATION. */

/*         LDT       INTEGER */
/*                   LDT IS THE LEADING DIMENSION OF THE ARRAY T. */

/*         N         INTEGER */
/*                   N IS THE ORDER OF THE SYSTEM. */

/*         B         DOUBLE PRECISION(N). */
/*                   B CONTAINS THE RIGHT HAND SIDE OF THE SYSTEM. */

/*         JOB       INTEGER */
/*                   JOB SPECIFIES WHAT KIND OF SYSTEM IS TO BE SOLVED. */
/*                   IF JOB IS */

/*                        00   SOLVE T*X=B, T LOWER TRIANGULAR, */
/*                        01   SOLVE T*X=B, T UPPER TRIANGULAR, */
/*                        10   SOLVE TRANS(T)*X=B, T LOWER TRIANGULAR, */
/*                        11   SOLVE TRANS(T)*X=B, T UPPER TRIANGULAR. */

/*     ON RETURN */

/*         B         B CONTAINS THE SOLUTION, IF INFO .EQ. 0. */
/*                   OTHERWISE B IS UNALTERED. */

/*         INFO      INTEGER */
/*                   INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR. */
/*                   OTHERWISE INFO CONTAINS THE INDEX OF */
/*                   THE FIRST ZERO DIAGONAL ELEMENT OF T. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G. W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SDOT */
/*     FORTRAN MOD */

/*     INTERNAL VARIABLES */


/*     BEGIN BLOCK PERMITTING ...EXITS TO 150 */

/*        CHECK FOR ZERO DIAGONAL ELEMENTS. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (*info = 1; *info <= i__1; ++(*info)) {
/*     ......EXIT */
	if (t[*info + *info * t_dim1] == 0.) {
	    goto L150;
	}
/* L10: */
    }
    *info = 0;

/*        DETERMINE THE TASK AND GO TO IT. */

    case__ = 1;
    if (*job % 10 != 0) {
	case__ = 2;
    }
    if (*job % 100 / 10 != 0) {
	case__ += 2;
    }
    switch (case__) {
	case 1:  goto L20;
	case 2:  goto L50;
	case 3:  goto L80;
	case 4:  goto L110;
    }

/*        SOLVE T*X=B FOR T LOWER TRIANGULAR */

L20:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L40;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	temp = -b[j - 1];
	i__2 = *n - j + 1;
	saxpy_(&i__2, &temp, &t[j + (j - 1) * t_dim1], &c__1, &b[j], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L30: */
    }
L40:
    goto L140;

/*        SOLVE T*X=B FOR T UPPER TRIANGULAR. */

L50:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L70;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	temp = -b[j + 1];
	saxpy_(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c__1, &b[1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L60: */
    }
L70:
    goto L140;

/*        SOLVE TRANS(T)*X=B FOR T LOWER TRIANGULAR. */

L80:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L100;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	i__2 = jj - 1;
	b[j] -= sdot_(&i__2, &t[j + 1 + j * t_dim1], &c__1, &b[j + 1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L90: */
    }
L100:
    goto L140;

/*        SOLVE TRANS(T)*X=B FOR T UPPER TRIANGULAR. */

L110:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L130;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	b[j] -= sdot_(&i__2, &t[j * t_dim1 + 1], &c__1, &b[1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L120: */
    }
L130:
L140:
L150:
    return 0;
} /* strsl_ */

/* Subroutine */ int strdi_(doublereal *t, integer *ldt, integer *n, 
	doublereal *det, integer *job, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, kb, km1, kp1;
    static doublereal ten, temp;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *), saxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);


/*     DTRDI COMPUTES THE DETERMINANT AND INVERSE OF A DOUBLE PRECISION */
/*     TRIANGULAR MATRIX. */

/*     ON ENTRY */

/*        T       DOUBLE PRECISION(LDT,N) */
/*                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO */
/*                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND */
/*                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE */
/*                USED TO STORE OTHER INFORMATION. */

/*        LDT     INTEGER */
/*                LDT IS THE LEADING DIMENSION OF THE ARRAY T. */

/*        N       INTEGER */
/*                N IS THE ORDER OF THE SYSTEM. */

/*        JOB     INTEGER */
/*                = 010       NO DET, INVERSE OF LOWER TRIANGULAR. */
/*                = 011       NO DET, INVERSE OF UPPER TRIANGULAR. */
/*                = 100       DET, NO INVERSE. */
/*                = 110       DET, INVERSE OF LOWER TRIANGULAR. */
/*                = 111       DET, INVERSE OF UPPER TRIANGULAR. */

/*     ON RETURN */

/*        T       INVERSE OF ORIGINAL MATRIX IF REQUESTED. */
/*                OTHERWISE UNCHANGED. */

/*        DET     DOUBLE PRECISION(2) */
/*                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED. */
/*                OTHERWISE NOT REFERENCED. */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. DABS(DET(1)) .LT. 10.0 */
/*                OR  DET(1) .EQ. 0.0 . */

/*        INFO    INTEGER */
/*                INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR */
/*                AND THE INVERSE IS REQUESTED. */
/*                OTHERWISE INFO CONTAINS THE INDEX OF */
/*                A ZERO DIAGONAL ELEMENT OF T. */


/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS SAXPY,SSCAL */
/*     FORTRAN DABS,MOD */

/*     INTERNAL VARIABLES */


/*     BEGIN BLOCK PERMITTING ...EXITS TO 180 */

/*        COMPUTE DETERMINANT */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --det;

    /* Function Body */
    if (*job / 100 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	det[1] = t[i__ + i__ * t_dim1] * det[1];
/*           ...EXIT */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (abs(det[1]) >= 1.) {
	    goto L20;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (abs(det[1]) < ten) {
	    goto L40;
	}
	det[1] /= ten;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*        COMPUTE INVERSE OF UPPER TRIANGULAR */

    if (*job / 10 % 10 == 0) {
	goto L170;
    }
    if (*job % 10 == 0) {
	goto L120;
    }
/*              BEGIN BLOCK PERMITTING ...EXITS TO 110 */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	*info = k;
/*              ......EXIT */
	if (t[k + k * t_dim1] == 0.) {
	    goto L110;
	}
	t[k + k * t_dim1] = 1. / t[k + k * t_dim1];
	temp = -t[k + k * t_dim1];
	i__2 = k - 1;
	sscal_(&i__2, &temp, &t[k * t_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    temp = t[k + j * t_dim1];
	    t[k + j * t_dim1] = 0.;
	    saxpy_(&k, &temp, &t[k * t_dim1 + 1], &c__1, &t[j * t_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }
    *info = 0;
L110:
    goto L160;
L120:

/*              COMPUTE INVERSE OF LOWER TRIANGULAR */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	*info = k;
/*     ............EXIT */
	if (t[k + k * t_dim1] == 0.) {
	    goto L180;
	}
	t[k + k * t_dim1] = 1. / t[k + k * t_dim1];
	temp = -t[k + k * t_dim1];
	if (k != *n) {
	    i__2 = *n - k;
	    sscal_(&i__2, &temp, &t[k + 1 + k * t_dim1], &c__1);
	}
	km1 = k - 1;
	if (km1 < 1) {
	    goto L140;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	    temp = t[k + j * t_dim1];
	    t[k + j * t_dim1] = 0.;
	    i__3 = *n - k + 1;
	    saxpy_(&i__3, &temp, &t[k + k * t_dim1], &c__1, &t[k + j * t_dim1]
		    , &c__1);
/* L130: */
	}
L140:
/* L150: */
	;
    }
    *info = 0;
L160:
L170:
L180:
    return 0;
} /* strdi_ */

/* Subroutine */ int sgtsl_(integer *n, doublereal *c__, doublereal *d__, 
	doublereal *e, doublereal *b, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer k;
    static doublereal t;
    static integer kb, kp1, nm1, nm2;


/*     DGTSL GIVEN A GENERAL TRIDIAGONAL MATRIX AND A RIGHT HAND */
/*     SIDE WILL FIND THE SOLUTION. */

/*     ON ENTRY */

/*        N       INTEGER */
/*                IS THE ORDER OF THE TRIDIAGONAL MATRIX. */

/*        C       DOUBLE PRECISION(N) */
/*                IS THE SUBDIAGONAL OF THE TRIDIAGONAL MATRIX. */
/*                C(2) THROUGH C(N) SHOULD CONTAIN THE SUBDIAGONAL. */
/*                ON OUTPUT C IS DESTROYED. */

/*        D       DOUBLE PRECISION(N) */
/*                IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX. */
/*                ON OUTPUT D IS DESTROYED. */

/*        E       DOUBLE PRECISION(N) */
/*                IS THE SUPERDIAGONAL OF THE TRIDIAGONAL MATRIX. */
/*                E(1) THROUGH E(N-1) SHOULD CONTAIN THE SUPERDIAGONAL. */
/*                ON OUTPUT E IS DESTROYED. */

/*        B       DOUBLE PRECISION(N) */
/*                IS THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       IS THE SOLUTION VECTOR. */

/*        INFO    INTEGER */
/*                = 0 NORMAL VALUE. */
/*                = K IF THE K-TH ELEMENT OF THE DIAGONAL BECOMES */
/*                    EXACTLY ZERO.  THE SUBROUTINE RETURNS WHEN */
/*                    THIS IS DETECTED. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JACK DONGARRA, ARGONNE NATIONAL LABORATORY. */

/*     NO EXTERNALS */
/*     FORTRAN DABS */

/*     INTERNAL VARIABLES */

/*     BEGIN BLOCK PERMITTING ...EXITS TO 100 */

    /* Parameter adjustments */
    --b;
    --e;
    --d__;
    --c__;

    /* Function Body */
    *info = 0;
    c__[1] = d__[1];
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L40;
    }
    d__[1] = e[1];
    e[1] = 0.;
    e[*n] = 0.;

    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*              FIND THE LARGEST OF THE TWO ROWS */

	if ((d__1 = c__[kp1], abs(d__1)) < (d__2 = c__[k], abs(d__2))) {
	    goto L10;
	}

/*                 INTERCHANGE ROW */

	t = c__[kp1];
	c__[kp1] = c__[k];
	c__[k] = t;
	t = d__[kp1];
	d__[kp1] = d__[k];
	d__[k] = t;
	t = e[kp1];
	e[kp1] = e[k];
	e[k] = t;
	t = b[kp1];
	b[kp1] = b[k];
	b[k] = t;
L10:

/*              ZERO ELEMENTS */

	if (c__[k] != 0.) {
	    goto L20;
	}
	*info = k;
/*     ............EXIT */
	goto L100;
L20:
	t = -c__[kp1] / c__[k];
	c__[kp1] = d__[kp1] + t * d__[k];
	d__[kp1] = e[kp1] + t * e[k];
	e[kp1] = 0.;
	b[kp1] += t * b[k];
/* L30: */
    }
L40:
    if (c__[*n] != 0.) {
	goto L50;
    }
    *info = *n;
    goto L90;
L50:

/*           BACK SOLVE */

    nm2 = *n - 2;
    b[*n] /= c__[*n];
    if (*n == 1) {
	goto L80;
    }
    b[nm1] = (b[nm1] - d__[nm1] * b[*n]) / c__[nm1];
    if (nm2 < 1) {
	goto L70;
    }
    i__1 = nm2;
    for (kb = 1; kb <= i__1; ++kb) {
	k = nm2 - kb + 1;
	b[k] = (b[k] - d__[k] * b[k + 1] - e[k] * b[k + 2]) / c__[k];
/* L60: */
    }
L70:
L80:
L90:
L100:

    return 0;
} /* sgtsl_ */

/* Subroutine */ int sptsl_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *b)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k;
    static doublereal t1, t2;
    static integer ke, kf, kp1, nm1, kbm1, nm1d2;


/*     DPTSL GIVEN A POSITIVE DEFINITE TRIDIAGONAL MATRIX AND A RIGHT */
/*     HAND SIDE WILL FIND THE SOLUTION. */

/*     ON ENTRY */

/*        N        INTEGER */
/*                 IS THE ORDER OF THE TRIDIAGONAL MATRIX. */

/*        D        DOUBLE PRECISION(N) */
/*                 IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX. */
/*                 ON OUTPUT D IS DESTROYED. */

/*        E        DOUBLE PRECISION(N) */
/*                 IS THE OFFDIAGONAL OF THE TRIDIAGONAL MATRIX. */
/*                 E(1) THROUGH E(N-1) SHOULD CONTAIN THE */
/*                 OFFDIAGONAL. */

/*        B        DOUBLE PRECISION(N) */
/*                 IS THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B        CONTAINS THE SOULTION. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JACK DONGARRA, ARGONNE NATIONAL LABORATORY. */

/*     NO EXTERNALS */
/*     FORTRAN MOD */

/*     INTERNAL VARIABLES */


/*     CHECK FOR 1 X 1 CASE */

    /* Parameter adjustments */
    --b;
    --e;
    --d__;

    /* Function Body */
    if (*n != 1) {
	goto L10;
    }
    b[1] /= d__[1];
    goto L70;
L10:
    nm1 = *n - 1;
    nm1d2 = nm1 / 2;
    if (*n == 2) {
	goto L30;
    }
    kbm1 = *n - 1;

/*           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF */
/*           SUPERDIAGONAL */

    i__1 = nm1d2;
    for (k = 1; k <= i__1; ++k) {
	t1 = e[k] / d__[k];
	d__[k + 1] -= t1 * e[k];
	b[k + 1] -= t1 * b[k];
	t2 = e[kbm1] / d__[kbm1 + 1];
	d__[kbm1] -= t2 * e[kbm1];
	b[kbm1] -= t2 * b[kbm1 + 1];
	--kbm1;
/* L20: */
    }
L30:
    kp1 = nm1d2 + 1;

/*        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER */

    if (*n % 2 != 0) {
	goto L40;
    }
    t1 = e[kp1] / d__[kp1];
    d__[kp1 + 1] -= t1 * e[kp1];
    b[kp1 + 1] -= t1 * b[kp1];
    ++kp1;
L40:

/*        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP */
/*        AND BOTTOM */

    b[kp1] /= d__[kp1];
    if (*n == 2) {
	goto L60;
    }
    k = kp1 - 1;
    ke = kp1 + nm1d2 - 1;
    i__1 = ke;
    for (kf = kp1; kf <= i__1; ++kf) {
	b[k] = (b[k] - e[k] * b[k + 1]) / d__[k];
	b[kf + 1] = (b[kf + 1] - e[kf] * b[kf]) / d__[kf + 1];
	--k;
/* L50: */
    }
L60:
    if (*n % 2 == 0) {
	b[1] = (b[1] - e[1] * b[2]) / d__[1];
    }
L70:
    return 0;
} /* sptsl_ */

/* Subroutine */ int schdc_(doublereal *a, integer *lda, integer *p, 
	doublereal *work, integer *jpvt, integer *job, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, k, l, kb, jp, pl, jt, pu, km1, kp1, plp1;
    static logical negk;
    static integer maxl;
    static doublereal temp;
    static logical swapk;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal maxdia;


/*     DCHDC COMPUTES THE CHOLESKY DECOMPOSITION OF A POSITIVE DEFINITE */
/*     MATRIX.  A PIVOTING OPTION ALLOWS THE USER TO ESTIMATE THE */
/*     CONDITION OF A POSITIVE DEFINITE MATRIX OR DETERMINE THE RANK */
/*     OF A POSITIVE SEMIDEFINITE MATRIX. */

/*     ON ENTRY */

/*         A      DOUBLE PRECISION(LDA,P). */
/*                A CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO */
/*                BE COMPUTED.  ONLT THE UPPER HALF OF A NEED BE STORED. */
/*                THE LOWER PART OF THE ARRAY A IS NOT REFERENCED. */

/*         LDA    INTEGER. */
/*                LDA IS THE LEADING DIMENSION OF THE ARRAY A. */

/*         P      INTEGER. */
/*                P IS THE ORDER OF THE MATRIX. */

/*         WORK   DOUBLE PRECISION. */
/*                WORK IS A WORK ARRAY. */

/*         JPVT   INTEGER(P). */
/*                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION */
/*                OF THE PIVOT ELEMENTS, IF PIVOTING HAS BEEN REQUESTED. */
/*                EACH DIAGONAL ELEMENT A(K,K) */
/*                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE */
/*                VALUE OF JPVT(K). */

/*                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL */
/*                                      ELEMENT. */

/*                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE ELEMENT. */

/*                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL ELEMENT. */

/*                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL ELEMENTS */
/*                ARE MOVED BY SYMMETRIC ROW AND COLUMN INTERCHANGES TO */
/*                THE BEGINNING OF THE ARRAY A AND FINAL */
/*                ELEMENTS TO THE END.  BOTH INITIAL AND FINAL ELEMENTS */
/*                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY */
/*                FREE ELEMENTS ARE MOVED.  AT THE K-TH STAGE OF THE */
/*                REDUCTION, IF A(K,K) IS OCCUPIED BY A FREE ELEMENT */
/*                IT IS INTERCHANGED WITH THE LARGEST FREE ELEMENT */
/*                A(L,L) WITH L .GE. K.  JPVT IS NOT REFERENCED IF */
/*                JOB .EQ. 0. */

/*        JOB     INTEGER. */
/*                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING. */
/*                IF JOB .EQ. 0, NO PIVOTING IS DONE. */
/*                IF JOB .NE. 0, PIVOTING IS DONE. */

/*     ON RETURN */

/*         A      A CONTAINS IN ITS UPPER HALF THE CHOLESKY FACTOR */
/*                OF THE MATRIX A AS IT HAS BEEN PERMUTED BY PIVOTING. */

/*         JPVT   JPVT(J) CONTAINS THE INDEX OF THE DIAGONAL ELEMENT */
/*                OF A THAT WAS MOVED INTO THE J-TH POSITION, */
/*                PROVIDED PIVOTING WAS REQUESTED. */

/*         INFO   CONTAINS THE INDEX OF THE LAST POSITIVE DIAGONAL */
/*                ELEMENT OF THE CHOLESKY FACTOR. */

/*     FOR POSITIVE DEFINITE MATRICES INFO = P IS THE NORMAL RETURN. */
/*     FOR PIVOTING WITH POSITIVE SEMIDEFINITE MATRICES INFO WILL */
/*     IN GENERAL BE LESS THAN P.  HOWEVER, INFO MAY BE GREATER THAN */
/*     THE RANK OF A, SINCE ROUNDING ERROR CAN CAUSE AN OTHERWISE ZERO */
/*     ELEMENT TO BE POSITIVE. INDEFINITE SYSTEMS WILL ALWAYS CAUSE */
/*     INFO TO BE LESS THAN P. */

/*     LINPACK. THIS VERSION DATED 03/19/79 . */
/*     J.J. DONGARRA AND G.W. STEWART, ARGONNE NATIONAL LABORATORY AND */
/*     UNIVERSITY OF MARYLAND. */


/*     BLAS SAXPY,SSWAP */
/*     FORTRAN DSQRT */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    --jpvt;

    /* Function Body */
    pl = 1;
    pu = 0;
    *info = *p;
    if (*job == 0) {
	goto L160;
    }

/*        PIVOTING HAS BEEN REQUESTED. REARRANGE THE */
/*        THE ELEMENTS ACCORDING TO JPVT. */

    i__1 = *p;
    for (k = 1; k <= i__1; ++k) {
	swapk = jpvt[k] > 0;
	negk = jpvt[k] < 0;
	jpvt[k] = k;
	if (negk) {
	    jpvt[k] = -jpvt[k];
	}
	if (! swapk) {
	    goto L60;
	}
	if (k == pl) {
	    goto L50;
	}
	i__2 = pl - 1;
	sswap_(&i__2, &a[k * a_dim1 + 1], &c__1, &a[pl * a_dim1 + 1], &c__1);
	temp = a[k + k * a_dim1];
	a[k + k * a_dim1] = a[pl + pl * a_dim1];
	a[pl + pl * a_dim1] = temp;
	plp1 = pl + 1;
	if (*p < plp1) {
	    goto L40;
	}
	i__2 = *p;
	for (j = plp1; j <= i__2; ++j) {
	    if (j >= k) {
		goto L10;
	    }
	    temp = a[pl + j * a_dim1];
	    a[pl + j * a_dim1] = a[j + k * a_dim1];
	    a[j + k * a_dim1] = temp;
	    goto L20;
L10:
	    if (j == k) {
		goto L20;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[pl + j * a_dim1];
	    a[pl + j * a_dim1] = temp;
L20:
/* L30: */
	    ;
	}
L40:
	jpvt[k] = jpvt[pl];
	jpvt[pl] = k;
L50:
	++pl;
L60:
/* L70: */
	;
    }
    pu = *p;
    if (*p < pl) {
	goto L150;
    }
    i__1 = *p;
    for (kb = pl; kb <= i__1; ++kb) {
	k = *p - kb + pl;
	if (jpvt[k] >= 0) {
	    goto L130;
	}
	jpvt[k] = -jpvt[k];
	if (pu == k) {
	    goto L120;
	}
	i__2 = k - 1;
	sswap_(&i__2, &a[k * a_dim1 + 1], &c__1, &a[pu * a_dim1 + 1], &c__1);
	temp = a[k + k * a_dim1];
	a[k + k * a_dim1] = a[pu + pu * a_dim1];
	a[pu + pu * a_dim1] = temp;
	kp1 = k + 1;
	if (*p < kp1) {
	    goto L110;
	}
	i__2 = *p;
	for (j = kp1; j <= i__2; ++j) {
	    if (j >= pu) {
		goto L80;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[j + pu * a_dim1];
	    a[j + pu * a_dim1] = temp;
	    goto L90;
L80:
	    if (j == pu) {
		goto L90;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[pu + j * a_dim1];
	    a[pu + j * a_dim1] = temp;
L90:
/* L100: */
	    ;
	}
L110:
	jt = jpvt[k];
	jpvt[k] = jpvt[pu];
	jpvt[pu] = jt;
L120:
	--pu;
L130:
/* L140: */
	;
    }
L150:
L160:
    i__1 = *p;
    for (k = 1; k <= i__1; ++k) {

/*        REDUCTION LOOP. */

	maxdia = a[k + k * a_dim1];
	kp1 = k + 1;
	maxl = k;

/*        DETERMINE THE PIVOT ELEMENT. */

	if (k < pl || k >= pu) {
	    goto L190;
	}
	i__2 = pu;
	for (l = kp1; l <= i__2; ++l) {
	    if (a[l + l * a_dim1] <= maxdia) {
		goto L170;
	    }
	    maxdia = a[l + l * a_dim1];
	    maxl = l;
L170:
/* L180: */
	    ;
	}
L190:

/*        QUIT IF THE PIVOT ELEMENT IS NOT POSITIVE. */

	if (maxdia > 0.) {
	    goto L200;
	}
	*info = k - 1;
/*     ......EXIT */
	goto L280;
L200:
	if (k == maxl) {
	    goto L210;
	}

/*           START THE PIVOTING AND UPDATE JPVT. */

	km1 = k - 1;
	sswap_(&km1, &a[k * a_dim1 + 1], &c__1, &a[maxl * a_dim1 + 1], &c__1);
	a[maxl + maxl * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = maxdia;
	jp = jpvt[maxl];
	jpvt[maxl] = jpvt[k];
	jpvt[k] = jp;
L210:

/*        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS. */

	work[k] = sqrt(a[k + k * a_dim1]);
	a[k + k * a_dim1] = work[k];
	if (*p < kp1) {
	    goto L260;
	}
	i__2 = *p;
	for (j = kp1; j <= i__2; ++j) {
	    if (k == maxl) {
		goto L240;
	    }
	    if (j >= maxl) {
		goto L220;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[j + maxl * a_dim1];
	    a[j + maxl * a_dim1] = temp;
	    goto L230;
L220:
	    if (j == maxl) {
		goto L230;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[maxl + j * a_dim1];
	    a[maxl + j * a_dim1] = temp;
L230:
L240:
	    a[k + j * a_dim1] /= work[k];
	    work[j] = a[k + j * a_dim1];
	    temp = -a[k + j * a_dim1];
	    i__3 = j - k;
	    saxpy_(&i__3, &temp, &work[kp1], &c__1, &a[kp1 + j * a_dim1], &
		    c__1);
/* L250: */
	}
L260:
/* L270: */
	;
    }
L280:
    return 0;
} /* schdc_ */

/* Subroutine */ int schud_(doublereal *r__, integer *ldr, integer *p, 
	doublereal *x, doublereal *z__, integer *ldz, integer *nz, doublereal 
	*y, doublereal *rho, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal t, xj;
    static integer jm1;
    static doublereal zeta, scale, azeta;
    extern /* Subroutine */ int srotg_(doublereal *, doublereal *, doublereal 
	    *, doublereal *);


/*     DCHUD UPDATES AN AUGMENTED CHOLESKY DECOMPOSITION OF THE */
/*     TRIANGULAR PART OF AN AUGMENTED QR DECOMPOSITION.  SPECIFICALLY, */
/*     GIVEN AN UPPER TRIANGULAR MATRIX R OF ORDER P, A ROW VECTOR */
/*     X, A COLUMN VECTOR Z, AND A SCALAR Y, DCHUD DETERMINES A */
/*     UNTIARY MATRIX U AND A SCALAR ZETA SUCH THAT */


/*                              (R  Z)     (RR   ZZ ) */
/*                         U  * (    )  =  (        ) , */
/*                              (X  Y)     ( 0  ZETA) */

/*     WHERE RR IS UPPER TRIANGULAR.  IF R AND Z HAVE BEEN */
/*     OBTAINED FROM THE FACTORIZATION OF A LEAST SQUARES */
/*     PROBLEM, THEN RR AND ZZ ARE THE FACTORS CORRESPONDING TO */
/*     THE PROBLEM WITH THE OBSERVATION (X,Y) APPENDED.  IN THIS */
/*     CASE, IF RHO IS THE NORM OF THE RESIDUAL VECTOR, THEN THE */
/*     NORM OF THE RESIDUAL VECTOR OF THE UPDATED PROBLEM IS */
/*     DSQRT(RHO**2 + ZETA**2).  DCHUD WILL SIMULTANEOUSLY UPDATE */
/*     SEVERAL TRIPLETS (Z,Y,RHO). */
/*     FOR A LESS TERSE DESCRIPTION OF WHAT DCHUD DOES AND HOW */
/*     IT MAY BE APPLIED, SEE THE LINPACK GUIDE. */

/*     THE MATRIX U IS DETERMINED AS THE PRODUCT U(P)*...*U(1), */
/*     WHERE U(I) IS A ROTATION IN THE (I,P+1) PLANE OF THE */
/*     FORM */

/*                       (     C(I)      S(I) ) */
/*                       (                    ) . */
/*                       (    -S(I)      C(I) ) */

/*     THE ROTATIONS ARE CHOSEN SO THAT C(I) IS DOUBLE PRECISION. */

/*     ON ENTRY */

/*         R      DOUBLE PRECISION(LDR,P), WHERE LDR .GE. P. */
/*                R CONTAINS THE UPPER TRIANGULAR MATRIX */
/*                THAT IS TO BE UPDATED.  THE PART OF R */
/*                BELOW THE DIAGONAL IS NOT REFERENCED. */

/*         LDR    INTEGER. */
/*                LDR IS THE LEADING DIMENSION OF THE ARRAY R. */

/*         P      INTEGER. */
/*                P IS THE ORDER OF THE MATRIX R. */

/*         X      DOUBLE PRECISION(P). */
/*                X CONTAINS THE ROW TO BE ADDED TO R.  X IS */
/*                NOT ALTERED BY DCHUD. */

/*         Z      DOUBLE PRECISION(LDZ,NZ), WHERE LDZ .GE. P. */
/*                Z IS AN ARRAY CONTAINING NZ P-VECTORS TO */
/*                BE UPDATED WITH R. */

/*         LDZ    INTEGER. */
/*                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z. */

/*         NZ     INTEGER. */
/*                NZ IS THE NUMBER OF VECTORS TO BE UPDATED */
/*                NZ MAY BE ZERO, IN WHICH CASE Z, Y, AND RHO */
/*                ARE NOT REFERENCED. */

/*         Y      DOUBLE PRECISION(NZ). */
/*                Y CONTAINS THE SCALARS FOR UPDATING THE VECTORS */
/*                Z.  Y IS NOT ALTERED BY DCHUD. */

/*         RHO    DOUBLE PRECISION(NZ). */
/*                RHO CONTAINS THE NORMS OF THE RESIDUAL */
/*                VECTORS THAT ARE TO BE UPDATED.  IF RHO(J) */
/*                IS NEGATIVE, IT IS LEFT UNALTERED. */

/*     ON RETURN */

/*         RC */
/*         RHO    CONTAIN THE UPDATED QUANTITIES. */
/*         Z */

/*         C      DOUBLE PRECISION(P). */
/*                C CONTAINS THE COSINES OF THE TRANSFORMING */
/*                ROTATIONS. */

/*         S      DOUBLE PRECISION(P). */
/*                S CONTAINS THE SINES OF THE TRANSFORMING */
/*                ROTATIONS. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     DCHUD USES THE FOLLOWING FUNCTIONS AND SUBROUTINES. */

/*     EXTENDED BLAS SROTG */
/*     FORTRAN DSQRT */


/*     UPDATE R. */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --x;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --y;
    --rho;
    --c__;
    --s;

    /* Function Body */
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	xj = x[j];

/*        APPLY THE PREVIOUS ROTATIONS. */

	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    t = c__[i__] * r__[i__ + j * r_dim1] + s[i__] * xj;
	    xj = c__[i__] * xj - s[i__] * r__[i__ + j * r_dim1];
	    r__[i__ + j * r_dim1] = t;
/* L10: */
	}
L20:

/*        COMPUTE THE NEXT ROTATION. */

	srotg_(&r__[j + j * r_dim1], &xj, &c__[j], &s[j]);
/* L30: */
    }

/*     IF REQUIRED, UPDATE Z AND RHO. */

    if (*nz < 1) {
	goto L70;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	zeta = y[j];
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    t = c__[i__] * z__[i__ + j * z_dim1] + s[i__] * zeta;
	    zeta = c__[i__] * zeta - s[i__] * z__[i__ + j * z_dim1];
	    z__[i__ + j * z_dim1] = t;
/* L40: */
	}
	azeta = abs(zeta);
	if (azeta == 0. || rho[j] < 0.) {
	    goto L50;
	}
	scale = azeta + rho[j];
/* Computing 2nd power */
	d__1 = azeta / scale;
/* Computing 2nd power */
	d__2 = rho[j] / scale;
	rho[j] = scale * sqrt(d__1 * d__1 + d__2 * d__2);
L50:
/* L60: */
	;
    }
L70:
    return 0;
} /* schud_ */

/* Subroutine */ int schdd_(doublereal *r__, integer *ldr, integer *p, 
	doublereal *x, doublereal *z__, integer *ldz, integer *nz, doublereal 
	*y, doublereal *rho, doublereal *c__, doublereal *s, integer *info)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal a, b;
    static integer i__, j;
    static doublereal t;
    static integer ii;
    static doublereal xx, zeta;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal norm;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static doublereal alpha;
    static real scale;
    static doublereal azeta;


/*     DCHDD DOWNDATES AN AUGMENTED CHOLESKY DECOMPOSITION OR THE */
/*     TRIANGULAR FACTOR OF AN AUGMENTED QR DECOMPOSITION. */
/*     SPECIFICALLY, GIVEN AN UPPER TRIANGULAR MATRIX R OF ORDER P,  A */
/*     ROW VECTOR X, A COLUMN VECTOR Z, AND A SCALAR Y, DCHDD */
/*     DETERMINEDS A ORTHOGONAL MATRIX U AND A SCALAR ZETA SUCH THAT */

/*                        (R   Z )     (RR  ZZ) */
/*                    U * (      )  =  (      ) , */
/*                        (0 ZETA)     ( X   Y) */

/*     WHERE RR IS UPPER TRIANGULAR.  IF R AND Z HAVE BEEN OBTAINED */
/*     FROM THE FACTORIZATION OF A LEAST SQUARES PROBLEM, THEN */
/*     RR AND ZZ ARE THE FACTORS CORRESPONDING TO THE PROBLEM */
/*     WITH THE OBSERVATION (X,Y) REMOVED.  IN THIS CASE, IF RHO */
/*     IS THE NORM OF THE RESIDUAL VECTOR, THEN THE NORM OF */
/*     THE RESIDUAL VECTOR OF THE DOWNDATED PROBLEM IS */
/*     DSQRT(RHO**2 - ZETA**2). DCHDD WILL SIMULTANEOUSLY DOWNDATE */
/*     SEVERAL TRIPLETS (Z,Y,RHO) ALONG WITH R. */
/*     FOR A LESS TERSE DESCRIPTION OF WHAT DCHDD DOES AND HOW */
/*     IT MAY BE APPLIED, SEE THE LINPACK GUIDE. */

/*     THE MATRIX U IS DETERMINED AS THE PRODUCT U(1)*...*U(P) */
/*     WHERE U(I) IS A ROTATION IN THE (P+1,I)-PLANE OF THE */
/*     FORM */

/*                       ( C(I)     -S(I)     ) */
/*                       (                    ) . */
/*                       ( S(I)       C(I)    ) */

/*     THE ROTATIONS ARE CHOSEN SO THAT C(I) IS DOUBLE PRECISION. */

/*     THE USER IS WARNED THAT A GIVEN DOWNDATING PROBLEM MAY */
/*     BE IMPOSSIBLE TO ACCOMPLISH OR MAY PRODUCE */
/*     INACCURATE RESULTS.  FOR EXAMPLE, THIS CAN HAPPEN */
/*     IF X IS NEAR A VECTOR WHOSE REMOVAL WILL REDUCE THE */
/*     RANK OF R.  BEWARE. */

/*     ON ENTRY */

/*         R      DOUBLE PRECISION(LDR,P), WHERE LDR .GE. P. */
/*                R CONTAINS THE UPPER TRIANGULAR MATRIX */
/*                THAT IS TO BE DOWNDATED.  THE PART OF  R */
/*                BELOW THE DIAGONAL IS NOT REFERENCED. */

/*         LDR    INTEGER. */
/*                LDR IS THE LEADING DIMENSION FO THE ARRAY R. */

/*         P      INTEGER. */
/*                P IS THE ORDER OF THE MATRIX R. */

/*         X      DOUBLE PRECISION(P). */
/*                X CONTAINS THE ROW VECTOR THAT IS TO */
/*                BE REMOVED FROM R.  X IS NOT ALTERED BY DCHDD. */

/*         Z      DOUBLE PRECISION(LDZ,NZ), WHERE LDZ .GE. P. */
/*                Z IS AN ARRAY OF NZ P-VECTORS WHICH */
/*                ARE TO BE DOWNDATED ALONG WITH R. */

/*         LDZ    INTEGER. */
/*                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z. */

/*         NZ     INTEGER. */
/*                NZ IS THE NUMBER OF VECTORS TO BE DOWNDATED */
/*                NZ MAY BE ZERO, IN WHICH CASE Z, Y, AND RHO */
/*                ARE NOT REFERENCED. */

/*         Y      DOUBLE PRECISION(NZ). */
/*                Y CONTAINS THE SCALARS FOR THE DOWNDATING */
/*                OF THE VECTORS Z.  Y IS NOT ALTERED BY DCHDD. */

/*         RHO    DOUBLE PRECISION(NZ). */
/*                RHO CONTAINS THE NORMS OF THE RESIDUAL */
/*                VECTORS THAT ARE TO BE DOWNDATED. */

/*     ON RETURN */

/*         R */
/*         Z      CONTAIN THE DOWNDATED QUANTITIES. */
/*         RHO */

/*         C      DOUBLE PRECISION(P). */
/*                C CONTAINS THE COSINES OF THE TRANSFORMING */
/*                ROTATIONS. */

/*         S      DOUBLE PRECISION(P). */
/*                S CONTAINS THE SINES OF THE TRANSFORMING */
/*                ROTATIONS. */

/*         INFO   INTEGER. */
/*                INFO IS SET AS FOLLOWS. */

/*                   INFO = 0  IF THE ENTIRE DOWNDATING */
/*                             WAS SUCCESSFUL. */

/*                   INFO =-1  IF R COULD NOT BE DOWNDATED. */
/*                             IN THIS CASE, ALL QUANTITIES */
/*                             ARE LEFT UNALTERED. */

/*                   INFO = 1  IF SOME RHO COULD NOT BE */
/*                             DOWNDATED.  THE OFFENDING RHOS ARE */
/*                             SET TO -1. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     DCHDD USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS. */

/*     FORTRAN DABS */
/*     BLAS SDOT, SNRM2 */


/*     SOLVE THE SYSTEM TRANS(R)*A = X, PLACING THE RESULT */
/*     IN THE ARRAY S. */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --x;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --y;
    --rho;
    --c__;
    --s;

    /* Function Body */
    *info = 0;
    s[1] = x[1] / r__[r_dim1 + 1];
    if (*p < 2) {
	goto L20;
    }
    i__1 = *p;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	s[j] = x[j] - sdot_(&i__2, &r__[j * r_dim1 + 1], &c__1, &s[1], &c__1);
	s[j] /= r__[j + j * r_dim1];
/* L10: */
    }
L20:
    norm = snrm2_(p, &s[1], &c__1);
    if (norm < 1.) {
	goto L30;
    }
    *info = -1;
    goto L120;
L30:
/* Computing 2nd power */
    d__1 = norm;
    alpha = sqrt(1. - d__1 * d__1);

/*        DETERMINE THE TRANSFORMATIONS. */

    i__1 = *p;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *p - ii + 1;
	scale = alpha + (d__1 = s[i__], abs(d__1));
	a = alpha / scale;
	b = s[i__] / scale;
/* Computing 2nd power */
	d__1 = a;
/* Computing 2nd power */
	d__2 = b;
	norm = sqrt(d__1 * d__1 + d__2 * d__2 + 0.);
	c__[i__] = a / norm;
	s[i__] = b / norm;
	alpha = scale * norm;
/* L40: */
    }

/*        APPLY THE TRANSFORMATIONS TO R. */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	xx = 0.;
	i__2 = j;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = j - ii + 1;
	    t = c__[i__] * xx + s[i__] * r__[i__ + j * r_dim1];
	    r__[i__ + j * r_dim1] = c__[i__] * r__[i__ + j * r_dim1] - s[i__] 
		    * xx;
	    xx = t;
/* L50: */
	}
/* L60: */
    }

/*        IF REQUIRED, DOWNDATE Z AND RHO. */

    if (*nz < 1) {
	goto L110;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	zeta = y[j];
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__ + j * z_dim1] = (z__[i__ + j * z_dim1] - s[i__] * zeta) / 
		    c__[i__];
	    zeta = c__[i__] * zeta - s[i__] * z__[i__ + j * z_dim1];
/* L70: */
	}
	azeta = abs(zeta);
	if (azeta <= rho[j]) {
	    goto L80;
	}
	*info = 1;
	rho[j] = -1.;
	goto L90;
L80:
/* Computing 2nd power */
	d__1 = azeta / rho[j];
	rho[j] *= sqrt(1. - d__1 * d__1);
L90:
/* L100: */
	;
    }
L110:
L120:
    return 0;
} /* schdd_ */

/* Subroutine */ int schex_(doublereal *r__, integer *ldr, integer *p, 
	integer *k, integer *l, doublereal *z__, integer *ldz, integer *nz, 
	doublereal *c__, doublereal *s, integer *job)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j;
    static doublereal t;
    static integer ii, jj, il, iu, km1, lm1, kp1, lmk;
    extern /* Subroutine */ int srotg_(doublereal *, doublereal *, doublereal 
	    *, doublereal *);


/*     DCHEX UPDATES THE CHOLESKY FACTORIZATION */

/*                   A = TRANS(R)*R */

/*     OF A POSITIVE DEFINITE MATRIX A OF ORDER P UNDER DIAGONAL */
/*     PERMUTATIONS OF THE FORM */

/*                   TRANS(E)*A*E */

/*     WHERE E IS A PERMUTATION MATRIX.  SPECIFICALLY, GIVEN */
/*     AN UPPER TRIANGULAR MATRIX R AND A PERMUTATION MATRIX */
/*     E (WHICH IS SPECIFIED BY K, L, AND JOB), DCHEX DETERMINES */
/*     A ORTHOGONAL MATRIX U SUCH THAT */

/*                           U*R*E = RR, */

/*     WHERE RR IS UPPER TRIANGULAR.  AT THE USERS OPTION, THE */
/*     TRANSFORMATION U WILL BE MULTIPLIED INTO THE ARRAY Z. */
/*     IF A = TRANS(X)*X, SO THAT R IS THE TRIANGULAR PART OF THE */
/*     QR FACTORIZATION OF X, THEN RR IS THE TRIANGULAR PART OF THE */
/*     QR FACTORIZATION OF X*E, I.E. X WITH ITS COLUMNS PERMUTED. */
/*     FOR A LESS TERSE DESCRIPTION OF WHAT DCHEX DOES AND HOW */
/*     IT MAY BE APPLIED, SEE THE LINPACK GUIDE. */

/*     THE MATRIX Q IS DETERMINED AS THE PRODUCT U(L-K)*...*U(1) */
/*     OF PLANE ROTATIONS OF THE FORM */

/*                           (    C(I)       S(I) ) */
/*                           (                    ) , */
/*                           (    -S(I)      C(I) ) */

/*     WHERE C(I) IS DOUBLE PRECISION, THE ROWS THESE ROTATIONS OPERATE */
/*     ON ARE DESCRIBED BELOW. */

/*     THERE ARE TWO TYPES OF PERMUTATIONS, WHICH ARE DETERMINED */
/*     BY THE VALUE OF JOB. */

/*     1. RIGHT CIRCULAR SHIFT (JOB = 1). */

/*         THE COLUMNS ARE REARRANGED IN THE FOLLOWING ORDER. */

/*                1,...,K-1,L,K,K+1,...,L-1,L+1,...,P. */

/*         U IS THE PRODUCT OF L-K ROTATIONS U(I), WHERE U(I) */
/*         ACTS IN THE (L-I,L-I+1)-PLANE. */

/*     2. LEFT CIRCULAR SHIFT (JOB = 2). */
/*         THE COLUMNS ARE REARRANGED IN THE FOLLOWING ORDER */

/*                1,...,K-1,K+1,K+2,...,L,K,L+1,...,P. */

/*         U IS THE PRODUCT OF L-K ROTATIONS U(I), WHERE U(I) */
/*         ACTS IN THE (K+I-1,K+I)-PLANE. */

/*     ON ENTRY */

/*         R      DOUBLE PRECISION(LDR,P), WHERE LDR.GE.P. */
/*                R CONTAINS THE UPPER TRIANGULAR FACTOR */
/*                THAT IS TO BE UPDATED.  ELEMENTS OF R */
/*                BELOW THE DIAGONAL ARE NOT REFERENCED. */

/*         LDR    INTEGER. */
/*                LDR IS THE LEADING DIMENSION OF THE ARRAY R. */

/*         P      INTEGER. */
/*                P IS THE ORDER OF THE MATRIX R. */

/*         K      INTEGER. */
/*                K IS THE FIRST COLUMN TO BE PERMUTED. */

/*         L      INTEGER. */
/*                L IS THE LAST COLUMN TO BE PERMUTED. */
/*                L MUST BE STRICTLY GREATER THAN K. */

/*         Z      DOUBLE PRECISION(LDZ,NZ), WHERE LDZ.GE.P. */
/*                Z IS AN ARRAY OF NZ P-VECTORS INTO WHICH THE */
/*                TRANSFORMATION U IS MULTIPLIED.  Z IS */
/*                NOT REFERENCED IF NZ = 0. */

/*         LDZ    INTEGER. */
/*                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z. */

/*         NZ     INTEGER. */
/*                NZ IS THE NUMBER OF COLUMNS OF THE MATRIX Z. */

/*         JOB    INTEGER. */
/*                JOB DETERMINES THE TYPE OF PERMUTATION. */
/*                       JOB = 1  RIGHT CIRCULAR SHIFT. */
/*                       JOB = 2  LEFT CIRCULAR SHIFT. */

/*     ON RETURN */

/*         R      CONTAINS THE UPDATED FACTOR. */

/*         Z      CONTAINS THE UPDATED MATRIX Z. */

/*         C      DOUBLE PRECISION(P). */
/*                C CONTAINS THE COSINES OF THE TRANSFORMING ROTATIONS. */

/*         S      DOUBLE PRECISION(P). */
/*                S CONTAINS THE SINES OF THE TRANSFORMING ROTATIONS. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     DCHEX USES THE FOLLOWING FUNCTIONS AND SUBROUTINES. */

/*     BLAS SROTG */
/*     FORTRAN MIN0 */


/*     INITIALIZE */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --c__;
    --s;

    /* Function Body */
    km1 = *k - 1;
    kp1 = *k + 1;
    lmk = *l - *k;
    lm1 = *l - 1;

/*     PERFORM THE APPROPRIATE TASK. */

    switch (*job) {
	case 1:  goto L10;
	case 2:  goto L130;
    }

/*     RIGHT CIRCULAR SHIFT. */

L10:

/*        REORDER THE COLUMNS. */

    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *l - i__ + 1;
	s[i__] = r__[ii + *l * r_dim1];
/* L20: */
    }
    i__1 = lm1;
    for (jj = *k; jj <= i__1; ++jj) {
	j = lm1 - jj + *k;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__ + (j + 1) * r_dim1] = r__[i__ + j * r_dim1];
/* L30: */
	}
	r__[j + 1 + (j + 1) * r_dim1] = 0.;
/* L40: */
    }
    if (*k == 1) {
	goto L60;
    }
    i__1 = km1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *l - i__ + 1;
	r__[i__ + *k * r_dim1] = s[ii];
/* L50: */
    }
L60:

/*        CALCULATE THE ROTATIONS. */

    t = s[1];
    i__1 = lmk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	srotg_(&s[i__ + 1], &t, &c__[i__], &s[i__]);
	t = s[i__ + 1];
/* L70: */
    }
    r__[*k + *k * r_dim1] = t;
    i__1 = *p;
    for (j = kp1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = 1, i__3 = *l - j + 1;
	il = max(i__2,i__3);
	i__2 = lmk;
	for (ii = il; ii <= i__2; ++ii) {
	    i__ = *l - ii;
	    t = c__[ii] * r__[i__ + j * r_dim1] + s[ii] * r__[i__ + 1 + j * 
		    r_dim1];
	    r__[i__ + 1 + j * r_dim1] = c__[ii] * r__[i__ + 1 + j * r_dim1] - 
		    s[ii] * r__[i__ + j * r_dim1];
	    r__[i__ + j * r_dim1] = t;
/* L80: */
	}
/* L90: */
    }

/*        IF REQUIRED, APPLY THE TRANSFORMATIONS TO Z. */

    if (*nz < 1) {
	goto L120;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = lmk;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = *l - ii;
	    t = c__[ii] * z__[i__ + j * z_dim1] + s[ii] * z__[i__ + 1 + j * 
		    z_dim1];
	    z__[i__ + 1 + j * z_dim1] = c__[ii] * z__[i__ + 1 + j * z_dim1] - 
		    s[ii] * z__[i__ + j * z_dim1];
	    z__[i__ + j * z_dim1] = t;
/* L100: */
	}
/* L110: */
    }
L120:
    goto L260;

/*     LEFT CIRCULAR SHIFT */

L130:

/*        REORDER THE COLUMNS */

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = lmk + i__;
	s[ii] = r__[i__ + *k * r_dim1];
/* L140: */
    }
    i__1 = lm1;
    for (j = *k; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__ + j * r_dim1] = r__[i__ + (j + 1) * r_dim1];
/* L150: */
	}
	jj = j - km1;
	s[jj] = r__[j + 1 + (j + 1) * r_dim1];
/* L160: */
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = lmk + i__;
	r__[i__ + *l * r_dim1] = s[ii];
/* L170: */
    }
    i__1 = *l;
    for (i__ = kp1; i__ <= i__1; ++i__) {
	r__[i__ + *l * r_dim1] = 0.;
/* L180: */
    }

/*        REDUCTION LOOP. */

    i__1 = *p;
    for (j = *k; j <= i__1; ++j) {
	if (j == *k) {
	    goto L200;
	}

/*              APPLY THE ROTATIONS. */

/* Computing MIN */
	i__2 = j - 1, i__3 = *l - 1;
	iu = min(i__2,i__3);
	i__2 = iu;
	for (i__ = *k; i__ <= i__2; ++i__) {
	    ii = i__ - *k + 1;
	    t = c__[ii] * r__[i__ + j * r_dim1] + s[ii] * r__[i__ + 1 + j * 
		    r_dim1];
	    r__[i__ + 1 + j * r_dim1] = c__[ii] * r__[i__ + 1 + j * r_dim1] - 
		    s[ii] * r__[i__ + j * r_dim1];
	    r__[i__ + j * r_dim1] = t;
/* L190: */
	}
L200:
	if (j >= *l) {
	    goto L210;
	}
	jj = j - *k + 1;
	t = s[jj];
	srotg_(&r__[j + j * r_dim1], &t, &c__[jj], &s[jj]);
L210:
/* L220: */
	;
    }

/*        APPLY THE ROTATIONS TO Z. */

    if (*nz < 1) {
	goto L250;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = lm1;
	for (i__ = *k; i__ <= i__2; ++i__) {
	    ii = i__ - km1;
	    t = c__[ii] * z__[i__ + j * z_dim1] + s[ii] * z__[i__ + 1 + j * 
		    z_dim1];
	    z__[i__ + 1 + j * z_dim1] = c__[ii] * z__[i__ + 1 + j * z_dim1] - 
		    s[ii] * z__[i__ + j * z_dim1];
	    z__[i__ + j * z_dim1] = t;
/* L230: */
	}
/* L240: */
    }
L250:
L260:
    return 0;
} /* schex_ */

/* Subroutine */ int sqrdc_(doublereal *x, integer *ldx, integer *n, integer *
	p, doublereal *qraux, integer *jpvt, doublereal *work, integer *job)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static integer j, l;
    static doublereal t;
    static integer jj, jp, pl, pu;
    static doublereal tt;
    static integer lp1, lup;
    static logical negj;
    static integer maxj;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *), snrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical swapj;
    extern /* Subroutine */ int sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal nrmxl;
    extern /* Subroutine */ int saxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal maxnrm;


/*     DQRDC USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR */
/*     FACTORIZATION OF AN N BY P MATRIX X.  COLUMN PIVOTING */
/*     BASED ON THE 2-NORMS OF THE REDUCED COLUMNS MAY BE */
/*     PERFORMED AT THE USERS OPTION. */

/*     ON ENTRY */

/*        X       DOUBLE PRECISION(LDX,P), WHERE LDX .GE. N. */
/*                X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE */
/*                COMPUTED. */

/*        LDX     INTEGER. */
/*                LDX IS THE LEADING DIMENSION OF THE ARRAY X. */

/*        N       INTEGER. */
/*                N IS THE NUMBER OF ROWS OF THE MATRIX X. */

/*        P       INTEGER. */
/*                P IS THE NUMBER OF COLUMNS OF THE MATRIX X. */

/*        JPVT    INTEGER(P). */
/*                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION */
/*                OF THE PIVOT COLUMNS.  THE K-TH COLUMN X(K) OF X */
/*                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE */
/*                VALUE OF JPVT(K). */

/*                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL */
/*                                      COLUMN. */

/*                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE COLUMN. */

/*                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL COLUMN. */

/*                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL COLUMNS */
/*                ARE MOVED TO THE BEGINNING OF THE ARRAY X AND FINAL */
/*                COLUMNS TO THE END.  BOTH INITIAL AND FINAL COLUMNS */
/*                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY */
/*                FREE COLUMNS ARE MOVED.  AT THE K-TH STAGE OF THE */
/*                REDUCTION, IF X(K) IS OCCUPIED BY A FREE COLUMN */
/*                IT IS INTERCHANGED WITH THE FREE COLUMN OF LARGEST */
/*                REDUCED NORM.  JPVT IS NOT REFERENCED IF */
/*                JOB .EQ. 0. */

/*        WORK    DOUBLE PRECISION(P). */
/*                WORK IS A WORK ARRAY.  WORK IS NOT REFERENCED IF */
/*                JOB .EQ. 0. */

/*        JOB     INTEGER. */
/*                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING. */
/*                IF JOB .EQ. 0, NO PIVOTING IS DONE. */
/*                IF JOB .NE. 0, PIVOTING IS DONE. */

/*     ON RETURN */

/*        X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER */
/*                TRIANGULAR MATRIX R OF THE QR FACTORIZATION. */
/*                BELOW ITS DIAGONAL X CONTAINS INFORMATION FROM */
/*                WHICH THE ORTHOGONAL PART OF THE DECOMPOSITION */
/*                CAN BE RECOVERED.  NOTE THAT IF PIVOTING HAS */
/*                BEEN REQUESTED, THE DECOMPOSITION IS NOT THAT */
/*                OF THE ORIGINAL MATRIX X BUT THAT OF X */
/*                WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT. */

/*        QRAUX   DOUBLE PRECISION(P). */
/*                QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER */
/*                THE ORTHOGONAL PART OF THE DECOMPOSITION. */

/*        JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE */
/*                ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO */
/*                THE K-TH COLUMN, IF PIVOTING WAS REQUESTED. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     DQRDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS. */

/*     BLAS SAXPY,SDOT,SSCAL,SSWAP,SNRM2 */
/*     FORTRAN DABS,DMAX1,MIN0,DSQRT */

/*     INTERNAL VARIABLES */



    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --qraux;
    --jpvt;
    --work;

    /* Function Body */
    pl = 1;
    pu = 0;
    if (*job == 0) {
	goto L60;
    }

/*        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS */
/*        ACCORDING TO JPVT. */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	swapj = jpvt[j] > 0;
	negj = jpvt[j] < 0;
	jpvt[j] = j;
	if (negj) {
	    jpvt[j] = -j;
	}
	if (! swapj) {
	    goto L10;
	}
	if (j != pl) {
	    sswap_(n, &x[pl * x_dim1 + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
	}
	jpvt[j] = jpvt[pl];
	jpvt[pl] = j;
	++pl;
L10:
/* L20: */
	;
    }
    pu = *p;
    i__1 = *p;
    for (jj = 1; jj <= i__1; ++jj) {
	j = *p - jj + 1;
	if (jpvt[j] >= 0) {
	    goto L40;
	}
	jpvt[j] = -jpvt[j];
	if (j == pu) {
	    goto L30;
	}
	sswap_(n, &x[pu * x_dim1 + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
	jp = jpvt[pu];
	jpvt[pu] = jpvt[j];
	jpvt[j] = jp;
L30:
	--pu;
L40:
/* L50: */
	;
    }
L60:

/*     COMPUTE THE NORMS OF THE FREE COLUMNS. */

    if (pu < pl) {
	goto L80;
    }
    i__1 = pu;
    for (j = pl; j <= i__1; ++j) {
	qraux[j] = snrm2_(n, &x[j * x_dim1 + 1], &c__1);
	work[j] = qraux[j];
/* L70: */
    }
L80:

/*     PERFORM THE HOUSEHOLDER REDUCTION OF X. */

    lup = min(*n,*p);
    i__1 = lup;
    for (l = 1; l <= i__1; ++l) {
	if (l < pl || l >= pu) {
	    goto L120;
	}

/*           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT */
/*           INTO THE PIVOT POSITION. */

	maxnrm = 0.;
	maxj = l;
	i__2 = pu;
	for (j = l; j <= i__2; ++j) {
	    if (qraux[j] <= maxnrm) {
		goto L90;
	    }
	    maxnrm = qraux[j];
	    maxj = j;
L90:
/* L100: */
	    ;
	}
	if (maxj == l) {
	    goto L110;
	}
	sswap_(n, &x[l * x_dim1 + 1], &c__1, &x[maxj * x_dim1 + 1], &c__1);
	qraux[maxj] = qraux[l];
	work[maxj] = work[l];
	jp = jpvt[maxj];
	jpvt[maxj] = jpvt[l];
	jpvt[l] = jp;
L110:
L120:
	qraux[l] = 0.;
	if (l == *n) {
	    goto L190;
	}

/*           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L. */

	i__2 = *n - l + 1;
	nrmxl = snrm2_(&i__2, &x[l + l * x_dim1], &c__1);
	if (nrmxl == 0.) {
	    goto L180;
	}
	if (x[l + l * x_dim1] != 0.) {
	    nrmxl = d_sign(&nrmxl, &x[l + l * x_dim1]);
	}
	i__2 = *n - l + 1;
	d__1 = 1. / nrmxl;
	sscal_(&i__2, &d__1, &x[l + l * x_dim1], &c__1);
	x[l + l * x_dim1] += 1.;

/*              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS, */
/*              UPDATING THE NORMS. */

	lp1 = l + 1;
	if (*p < lp1) {
	    goto L170;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    t = -sdot_(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1) / x[l + l * x_dim1];
	    i__3 = *n - l + 1;
	    saxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
	    if (j < pl || j > pu) {
		goto L150;
	    }
	    if (qraux[j] == 0.) {
		goto L150;
	    }
/* Computing 2nd power */
	    d__2 = (d__1 = x[l + j * x_dim1], abs(d__1)) / qraux[j];
	    tt = 1. - d__2 * d__2;
	    tt = max(tt,0.);
	    t = tt;
/* Computing 2nd power */
	    d__1 = qraux[j] / work[j];
	    tt = tt * .05 * (d__1 * d__1) + 1.;
	    if (tt == 1.) {
		goto L130;
	    }
	    qraux[j] *= sqrt(t);
	    goto L140;
L130:
	    i__3 = *n - l;
	    qraux[j] = snrm2_(&i__3, &x[l + 1 + j * x_dim1], &c__1);
	    work[j] = qraux[j];
L140:
L150:
/* L160: */
	    ;
	}
L170:

/*              SAVE THE TRANSFORMATION. */

	qraux[l] = x[l + l * x_dim1];
	x[l + l * x_dim1] = -nrmxl;
L180:
L190:
/* L200: */
	;
    }
    return 0;
} /* sqrdc_ */

/* Subroutine */ int sqrsl_(doublereal *x, integer *ldx, integer *n, integer *
	k, doublereal *qraux, doublereal *y, doublereal *qy, doublereal *qty, 
	doublereal *b, doublereal *rsd, doublereal *xb, integer *job, integer 
	*info)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    static doublereal t;
    static logical cb;
    static integer jj;
    static logical cr;
    static integer ju, kp1;
    static logical cxb, cqy;
    static doublereal temp;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical cqty;
    extern /* Subroutine */ int scopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);


/*     DQRSL APPLIES THE OUTPUT OF DQRDC TO COMPUTE COORDINATE */
/*     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS. */
/*     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX */

/*            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K))) */

/*     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL */
/*     N X P MATRIX X THAT WAS INPUT TO DQRDC (IF NO PIVOTING WAS */
/*     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR */
/*     ORIGINAL ORDER).  DQRDC PRODUCES A FACTORED ORTHOGONAL MATRIX Q */
/*     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT */

/*              XK = Q * (R) */
/*                       (0) */

/*     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS */
/*     X AND QRAUX. */

/*     ON ENTRY */

/*        X      DOUBLE PRECISION(LDX,P). */
/*               X CONTAINS THE OUTPUT OF DQRDC. */

/*        LDX    INTEGER. */
/*               LDX IS THE LEADING DIMENSION OF THE ARRAY X. */

/*        N      INTEGER. */
/*               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST */
/*               HAVE THE SAME VALUE AS N IN DQRDC. */

/*        K      INTEGER. */
/*               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K */
/*               MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE */
/*               SAME AS IN THE CALLING SEQUENCE TO DQRDC. */

/*        QRAUX  DOUBLE PRECISION(P). */
/*               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM DQRDC. */

/*        Y      DOUBLE PRECISION(N) */
/*               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED */
/*               BY DQRSL. */

/*        JOB    INTEGER. */
/*               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS */
/*               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING */
/*               MEANING. */

/*                    IF A.NE.0, COMPUTE QY. */
/*                    IF B,C,D, OR E .NE. 0, COMPUTE QTY. */
/*                    IF C.NE.0, COMPUTE B. */
/*                    IF D.NE.0, COMPUTE RSD. */
/*                    IF E.NE.0, COMPUTE XB. */

/*               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB */
/*               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR */
/*               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING */
/*               SEQUENCE. */

/*     ON RETURN */

/*        QY     DOUBLE PRECISION(N). */
/*               QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN */
/*               REQUESTED. */

/*        QTY    DOUBLE PRECISION(N). */
/*               QTY CONTAINS TRANS(Q)*Y, IF ITS COMPUTATION HAS */
/*               BEEN REQUESTED.  HERE TRANS(Q) IS THE */
/*               TRANSPOSE OF THE MATRIX Q. */

/*        B      DOUBLE PRECISION(K) */
/*               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM */

/*                    MINIMIZE NORM2(Y - XK*B), */

/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT */
/*               IF PIVOTING WAS REQUESTED IN DQRDC, THE J-TH */
/*               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J) */
/*               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO DQRDC.) */

/*        RSD    DOUBLE PRECISION(N). */
/*               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B, */
/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS */
/*               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE */
/*               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK. */

/*        XB     DOUBLE PRECISION(N). */
/*               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B, */
/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO */
/*               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE */
/*               OF X. */

/*        INFO   INTEGER. */
/*               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS */
/*               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN */
/*               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO */
/*               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED. */

/*     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED */
/*     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE */
/*     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM. */
/*     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME */
/*     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A */
/*     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE */
/*     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS */
/*     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE */
/*     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE */
/*     COMPUTED.  THUS THE CALLING SEQUENCE */

/*          CALL DQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO) */

/*     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD */
/*     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING */
/*     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR */
/*     A SINGLE CALLINNG SEQUENCE. */

/*          1. (Y,QTY,B) (RSD) (XB) (QY) */

/*          2. (Y,QTY,RSD) (B) (XB) (QY) */

/*          3. (Y,QTY,XB) (B) (RSD) (QY) */

/*          4. (Y,QY) (QTY,B) (RSD) (XB) */

/*          5. (Y,QY) (QTY,RSD) (B) (XB) */

/*          6. (Y,QY) (QTY,XB) (B) (RSD) */

/*     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO */
/*     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     DQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS. */

/*     BLAS SAXPY,SCOPY,SDOT */
/*     FORTRAN DABS,MIN0,MOD */

/*     INTERNAL VARIABLES */



/*     SET INFO FLAG. */

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --qraux;
    --y;
    --qy;
    --qty;
    --b;
    --rsd;
    --xb;

    /* Function Body */
    *info = 0;

/*     DETERMINE WHAT IS TO BE COMPUTED. */

    cqy = *job / 10000 != 0;
    cqty = *job % 10000 != 0;
    cb = *job % 1000 / 100 != 0;
    cr = *job % 100 / 10 != 0;
    cxb = *job % 10 != 0;
/* Computing MIN */
    i__1 = *k, i__2 = *n - 1;
    ju = min(i__1,i__2);

/*     SPECIAL ACTION WHEN N=1. */

    if (ju != 0) {
	goto L40;
    }
    if (cqy) {
	qy[1] = y[1];
    }
    if (cqty) {
	qty[1] = y[1];
    }
    if (cxb) {
	xb[1] = y[1];
    }
    if (! cb) {
	goto L30;
    }
    if (x[x_dim1 + 1] != 0.) {
	goto L10;
    }
    *info = 1;
    goto L20;
L10:
    b[1] = y[1] / x[x_dim1 + 1];
L20:
L30:
    if (cr) {
	rsd[1] = 0.;
    }
    goto L250;
L40:

/*        SET UP TO COMPUTE QY OR QTY. */

    if (cqy) {
	scopy_(n, &y[1], &c__1, &qy[1], &c__1);
    }
    if (cqty) {
	scopy_(n, &y[1], &c__1, &qty[1], &c__1);
    }
    if (! cqy) {
	goto L70;
    }

/*           COMPUTE QY. */

    i__1 = ju;
    for (jj = 1; jj <= i__1; ++jj) {
	j = ju - jj + 1;
	if (qraux[j] == 0.) {
	    goto L50;
	}
	temp = x[j + j * x_dim1];
	x[j + j * x_dim1] = qraux[j];
	i__2 = *n - j + 1;
	t = -sdot_(&i__2, &x[j + j * x_dim1], &c__1, &qy[j], &c__1) / x[j + j 
		* x_dim1];
	i__2 = *n - j + 1;
	saxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qy[j], &c__1);
	x[j + j * x_dim1] = temp;
L50:
/* L60: */
	;
    }
L70:
    if (! cqty) {
	goto L100;
    }

/*           COMPUTE TRANS(Q)*Y. */

    i__1 = ju;
    for (j = 1; j <= i__1; ++j) {
	if (qraux[j] == 0.) {
	    goto L80;
	}
	temp = x[j + j * x_dim1];
	x[j + j * x_dim1] = qraux[j];
	i__2 = *n - j + 1;
	t = -sdot_(&i__2, &x[j + j * x_dim1], &c__1, &qty[j], &c__1) / x[j + 
		j * x_dim1];
	i__2 = *n - j + 1;
	saxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qty[j], &c__1);
	x[j + j * x_dim1] = temp;
L80:
/* L90: */
	;
    }
L100:

/*        SET UP TO COMPUTE B, RSD, OR XB. */

    if (cb) {
	scopy_(k, &qty[1], &c__1, &b[1], &c__1);
    }
    kp1 = *k + 1;
    if (cxb) {
	scopy_(k, &qty[1], &c__1, &xb[1], &c__1);
    }
    if (cr && *k < *n) {
	i__1 = *n - *k;
	scopy_(&i__1, &qty[kp1], &c__1, &rsd[kp1], &c__1);
    }
    if (! cxb || kp1 > *n) {
	goto L120;
    }
    i__1 = *n;
    for (i__ = kp1; i__ <= i__1; ++i__) {
	xb[i__] = 0.;
/* L110: */
    }
L120:
    if (! cr) {
	goto L140;
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rsd[i__] = 0.;
/* L130: */
    }
L140:
    if (! cb) {
	goto L190;
    }

/*           COMPUTE B. */

    i__1 = *k;
    for (jj = 1; jj <= i__1; ++jj) {
	j = *k - jj + 1;
	if (x[j + j * x_dim1] != 0.) {
	    goto L150;
	}
	*info = j;
/*           ......EXIT */
	goto L180;
L150:
	b[j] /= x[j + j * x_dim1];
	if (j == 1) {
	    goto L160;
	}
	t = -b[j];
	i__2 = j - 1;
	saxpy_(&i__2, &t, &x[j * x_dim1 + 1], &c__1, &b[1], &c__1);
L160:
/* L170: */
	;
    }
L180:
L190:
    if (! cr && ! cxb) {
	goto L240;
    }

/*           COMPUTE RSD OR XB AS REQUIRED. */

    i__1 = ju;
    for (jj = 1; jj <= i__1; ++jj) {
	j = ju - jj + 1;
	if (qraux[j] == 0.) {
	    goto L220;
	}
	temp = x[j + j * x_dim1];
	x[j + j * x_dim1] = qraux[j];
	if (! cr) {
	    goto L200;
	}
	i__2 = *n - j + 1;
	t = -sdot_(&i__2, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1) / x[j + 
		j * x_dim1];
	i__2 = *n - j + 1;
	saxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1);
L200:
	if (! cxb) {
	    goto L210;
	}
	i__2 = *n - j + 1;
	t = -sdot_(&i__2, &x[j + j * x_dim1], &c__1, &xb[j], &c__1) / x[j + j 
		* x_dim1];
	i__2 = *n - j + 1;
	saxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &xb[j], &c__1);
L210:
	x[j + j * x_dim1] = temp;
L220:
/* L230: */
	;
    }
L240:
L250:
    return 0;
} /* sqrsl_ */

/* Subroutine */ int ssvdc_(doublereal *x, integer *ldx, integer *n, integer *
	p, doublereal *s, doublereal *e, doublereal *u, integer *ldu, 
	doublereal *v, integer *ldv, doublereal *work, integer *job, integer *
	info)
{
    /* System generated locals */
    integer x_dim1, x_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal b, c__, f, g;
    static integer i__, j, k, l, m;
    static doublereal t, t1, el;
    static integer kk;
    static doublereal cs;
    static integer ll, mm, ls;
    static doublereal sl;
    static integer lu;
    static doublereal sm, sn;
    static integer lm1, mm1, lp1, mp1, nct, ncu, lls, nrt;
    static doublereal emm1, smm1;
    static integer kase, jobu, iter;
    extern doublereal sdot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal test;
    extern /* Subroutine */ int srot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer nctp1;
    extern doublereal snrm2_(integer *, doublereal *, integer *);
    static integer nrtp1;
    static doublereal scale;
    extern /* Subroutine */ int sscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal shift;
    static integer maxit;
    static logical wantu, wantv;
    extern /* Subroutine */ int srotg_(doublereal *, doublereal *, doublereal 
	    *, doublereal *), sswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), saxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal ztest;



/*     DSVDC IS A SUBROUTINE TO REDUCE A DOUBLE PRECISION NXP MATRIX X */
/*     BY ORTHOGONAL TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE */
/*     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE */
/*     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS, */
/*     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS. */

/*     ON ENTRY */

/*         X         DOUBLE PRECISION(LDX,P), WHERE LDX.GE.N. */
/*                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE */
/*                   DECOMPOSITION IS TO BE COMPUTED.  X IS */
/*                   DESTROYED BY DSVDC. */

/*         LDX       INTEGER. */
/*                   LDX IS THE LEADING DIMENSION OF THE ARRAY X. */

/*         N         INTEGER. */
/*                   N IS THE NUMBER OF COLUMNS OF THE MATRIX X. */

/*         P         INTEGER. */
/*                   P IS THE NUMBER OF ROWS OF THE MATRIX X. */

/*         LDU       INTEGER. */
/*                   LDU IS THE LEADING DIMENSION OF THE ARRAY U. */
/*                   (SEE BELOW). */

/*         LDV       INTEGER. */
/*                   LDV IS THE LEADING DIMENSION OF THE ARRAY V. */
/*                   (SEE BELOW). */

/*         WORK      DOUBLE PRECISION(N). */
/*                   WORK IS A SCRATCH ARRAY. */

/*         JOB       INTEGER. */
/*                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR */
/*                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB */
/*                   WITH THE FOLLOWING MEANING */

/*                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR */
/*                                  VECTORS. */
/*                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS */
/*                                  IN U. */
/*                        A.GE.2    RETURN THE FIRST MIN(N,P) SINGULAR */
/*                                  VECTORS IN U. */
/*                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR */
/*                                  VECTORS. */
/*                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS */
/*                                  IN V. */

/*     ON RETURN */

/*         S         DOUBLE PRECISION(MM), WHERE MM=MIN(N+1,P). */
/*                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE */
/*                   SINGULAR VALUES OF X ARRANGED IN DESCENDING */
/*                   ORDER OF MAGNITUDE. */

/*         E         DOUBLE PRECISION(P). */
/*                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE */
/*                   DISCUSSION OF INFO FOR EXCEPTIONS. */

/*         U         DOUBLE PRECISION(LDU,K), WHERE LDU.GE.N.  IF */
/*                                   JOBA.EQ.1 THEN K.EQ.N, IF JOBA.GE.2 */
/*                                   THEN K.EQ.MIN(N,P). */
/*                   U CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS. */
/*                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P */
/*                   OR IF JOBA.EQ.2, THEN U MAY BE IDENTIFIED WITH X */
/*                   IN THE SUBROUTINE CALL. */

/*         V         DOUBLE PRECISION(LDV,P), WHERE LDV.GE.P. */
/*                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS. */
/*                   V IS NOT REFERENCED IF JOB.EQ.0.  IF P.LE.N, */
/*                   THEN V MAY BE IDENTIFIED WITH X IN THE */
/*                   SUBROUTINE CALL. */

/*         INFO      INTEGER. */
/*                   THE SINGULAR VALUES (AND THEIR CORRESPONDING */
/*                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M) */
/*                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF */
/*                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR */
/*                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX */
/*                   B = TRANS(U)*X*V IS THE BIDIAGONAL MATRIX */
/*                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE */
/*                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (TRANS(U) */
/*                   IS THE TRANSPOSE OF U).  THUS THE SINGULAR */
/*                   VALUES OF X AND B ARE THE SAME. */

/*     LINPACK. THIS VERSION DATED 03/19/79 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     DSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS. */

/*     EXTERNAL SROT */
/*     BLAS SAXPY,SDOT,SSCAL,SSWAP,SNRM2,SROTG */
/*     FORTRAN DABS,DMAX1,MAX0,MIN0,MOD,DSQRT */

/*     INTERNAL VARIABLES */



/*     SET THE MAXIMUM NUMBER OF ITERATIONS. */

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --s;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --work;

    /* Function Body */
    maxit = 30;

/*     DETERMINE WHAT IS TO BE COMPUTED. */

    wantu = FALSE_;
    wantv = FALSE_;
    jobu = *job % 100 / 10;
    ncu = *n;
    if (jobu > 1) {
	ncu = min(*n,*p);
    }
    if (jobu != 0) {
	wantu = TRUE_;
    }
    if (*job % 10 != 0) {
	wantv = TRUE_;
    }

/*     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS */
/*     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E. */

    *info = 0;
/* Computing MIN */
    i__1 = *n - 1;
    nct = min(i__1,*p);
/* Computing MAX */
/* Computing MIN */
    i__3 = *p - 2;
    i__1 = 0, i__2 = min(i__3,*n);
    nrt = max(i__1,i__2);
    lu = max(nct,nrt);
    if (lu < 1) {
	goto L170;
    }
    i__1 = lu;
    for (l = 1; l <= i__1; ++l) {
	lp1 = l + 1;
	if (l > nct) {
	    goto L20;
	}

/*           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND */
/*           PLACE THE L-TH DIAGONAL IN S(L). */

	i__2 = *n - l + 1;
	s[l] = snrm2_(&i__2, &x[l + l * x_dim1], &c__1);
	if (s[l] == 0.) {
	    goto L10;
	}
	if (x[l + l * x_dim1] != 0.) {
	    s[l] = d_sign(&s[l], &x[l + l * x_dim1]);
	}
	i__2 = *n - l + 1;
	d__1 = 1. / s[l];
	sscal_(&i__2, &d__1, &x[l + l * x_dim1], &c__1);
	x[l + l * x_dim1] += 1.;
L10:
	s[l] = -s[l];
L20:
	if (*p < lp1) {
	    goto L50;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    if (l > nct) {
		goto L30;
	    }
	    if (s[l] == 0.) {
		goto L30;
	    }

/*              APPLY THE TRANSFORMATION. */

	    i__3 = *n - l + 1;
	    t = -sdot_(&i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1) / x[l + l * x_dim1];
	    i__3 = *n - l + 1;
	    saxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
L30:

/*           PLACE THE L-TH ROW OF X INTO  E FOR THE */
/*           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION. */

	    e[j] = x[l + j * x_dim1];
/* L40: */
	}
L50:
	if (! wantu || l > nct) {
	    goto L70;
	}

/*           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK */
/*           MULTIPLICATION. */

	i__2 = *n;
	for (i__ = l; i__ <= i__2; ++i__) {
	    u[i__ + l * u_dim1] = x[i__ + l * x_dim1];
/* L60: */
	}
L70:
	if (l > nrt) {
	    goto L150;
	}

/*           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE */
/*           L-TH SUPER-DIAGONAL IN E(L). */

	i__2 = *p - l;
	e[l] = snrm2_(&i__2, &e[lp1], &c__1);
	if (e[l] == 0.) {
	    goto L80;
	}
	if (e[lp1] != 0.) {
	    e[l] = d_sign(&e[l], &e[lp1]);
	}
	i__2 = *p - l;
	d__1 = 1. / e[l];
	sscal_(&i__2, &d__1, &e[lp1], &c__1);
	e[lp1] += 1.;
L80:
	e[l] = -e[l];
	if (lp1 > *n || e[l] == 0.) {
	    goto L120;
	}

/*              APPLY THE TRANSFORMATION. */

	i__2 = *n;
	for (i__ = lp1; i__ <= i__2; ++i__) {
	    work[i__] = 0.;
/* L90: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    saxpy_(&i__3, &e[j], &x[lp1 + j * x_dim1], &c__1, &work[lp1], &
		    c__1);
/* L100: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    d__1 = -e[j] / e[lp1];
	    saxpy_(&i__3, &d__1, &work[lp1], &c__1, &x[lp1 + j * x_dim1], &
		    c__1);
/* L110: */
	}
L120:
	if (! wantv) {
	    goto L140;
	}

/*              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT */
/*              BACK MULTIPLICATION. */

	i__2 = *p;
	for (i__ = lp1; i__ <= i__2; ++i__) {
	    v[i__ + l * v_dim1] = e[i__];
/* L130: */
	}
L140:
L150:
/* L160: */
	;
    }
L170:

/*     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M. */

/* Computing MIN */
    i__1 = *p, i__2 = *n + 1;
    m = min(i__1,i__2);
    nctp1 = nct + 1;
    nrtp1 = nrt + 1;
    if (nct < *p) {
	s[nctp1] = x[nctp1 + nctp1 * x_dim1];
    }
    if (*n < m) {
	s[m] = 0.;
    }
    if (nrtp1 < m) {
	e[nrtp1] = x[nrtp1 + m * x_dim1];
    }
    e[m] = 0.;

/*     IF REQUIRED, GENERATE U. */

    if (! wantu) {
	goto L300;
    }
    if (ncu < nctp1) {
	goto L200;
    }
    i__1 = ncu;
    for (j = nctp1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    u[i__ + j * u_dim1] = 0.;
/* L180: */
	}
	u[j + j * u_dim1] = 1.;
/* L190: */
    }
L200:
    if (nct < 1) {
	goto L290;
    }
    i__1 = nct;
    for (ll = 1; ll <= i__1; ++ll) {
	l = nct - ll + 1;
	if (s[l] == 0.) {
	    goto L250;
	}
	lp1 = l + 1;
	if (ncu < lp1) {
	    goto L220;
	}
	i__2 = ncu;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    t = -sdot_(&i__3, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &
		    c__1) / u[l + l * u_dim1];
	    i__3 = *n - l + 1;
	    saxpy_(&i__3, &t, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &
		    c__1);
/* L210: */
	}
L220:
	i__2 = *n - l + 1;
	sscal_(&i__2, &c_b1033, &u[l + l * u_dim1], &c__1);
	u[l + l * u_dim1] += 1.;
	lm1 = l - 1;
	if (lm1 < 1) {
	    goto L240;
	}
	i__2 = lm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    u[i__ + l * u_dim1] = 0.;
/* L230: */
	}
L240:
	goto L270;
L250:
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    u[i__ + l * u_dim1] = 0.;
/* L260: */
	}
	u[l + l * u_dim1] = 1.;
L270:
/* L280: */
	;
    }
L290:
L300:

/*     IF IT IS REQUIRED, GENERATE V. */

    if (! wantv) {
	goto L350;
    }
    i__1 = *p;
    for (ll = 1; ll <= i__1; ++ll) {
	l = *p - ll + 1;
	lp1 = l + 1;
	if (l > nrt) {
	    goto L320;
	}
	if (e[l] == 0.) {
	    goto L320;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *p - l;
	    t = -sdot_(&i__3, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1) / v[lp1 + l * v_dim1];
	    i__3 = *p - l;
	    saxpy_(&i__3, &t, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1);
/* L310: */
	}
L320:
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    v[i__ + l * v_dim1] = 0.;
/* L330: */
	}
	v[l + l * v_dim1] = 1.;
/* L340: */
    }
L350:

/*     MAIN ITERATION LOOP FOR THE SINGULAR VALUES. */

    mm = m;
    iter = 0;
L360:

/*        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND. */

/*     ...EXIT */
    if (m == 0) {
	goto L620;
    }

/*        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET */
/*        FLAG AND RETURN. */

    if (iter < maxit) {
	goto L370;
    }
    *info = m;
/*     ......EXIT */
    goto L620;
L370:

/*        THIS SECTION OF THE PROGRAM INSPECTS FOR */
/*        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON */
/*        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS. */

/*           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M */
/*           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M */
/*           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND */
/*                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP). */
/*           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE). */

    i__1 = m;
    for (ll = 1; ll <= i__1; ++ll) {
	l = m - ll;
/*        ...EXIT */
	if (l == 0) {
	    goto L400;
	}
	test = (d__1 = s[l], abs(d__1)) + (d__2 = s[l + 1], abs(d__2));
	ztest = test + (d__1 = e[l], abs(d__1));
	if (ztest != test) {
	    goto L380;
	}
	e[l] = 0.;
/*        ......EXIT */
	goto L400;
L380:
/* L390: */
	;
    }
L400:
    if (l != m - 1) {
	goto L410;
    }
    kase = 4;
    goto L480;
L410:
    lp1 = l + 1;
    mp1 = m + 1;
    i__1 = mp1;
    for (lls = lp1; lls <= i__1; ++lls) {
	ls = m - lls + lp1;
/*           ...EXIT */
	if (ls == l) {
	    goto L440;
	}
	test = 0.;
	if (ls != m) {
	    test += (d__1 = e[ls], abs(d__1));
	}
	if (ls != l + 1) {
	    test += (d__1 = e[ls - 1], abs(d__1));
	}
	ztest = test + (d__1 = s[ls], abs(d__1));
	if (ztest != test) {
	    goto L420;
	}
	s[ls] = 0.;
/*           ......EXIT */
	goto L440;
L420:
/* L430: */
	;
    }
L440:
    if (ls != l) {
	goto L450;
    }
    kase = 3;
    goto L470;
L450:
    if (ls != m) {
	goto L460;
    }
    kase = 1;
    goto L470;
L460:
    kase = 2;
    l = ls;
L470:
L480:
    ++l;

/*        PERFORM THE TASK INDICATED BY KASE. */

    switch (kase) {
	case 1:  goto L490;
	case 2:  goto L520;
	case 3:  goto L540;
	case 4:  goto L570;
    }

/*        DEFLATE NEGLIGIBLE S(M). */

L490:
    mm1 = m - 1;
    f = e[m - 1];
    e[m - 1] = 0.;
    i__1 = mm1;
    for (kk = l; kk <= i__1; ++kk) {
	k = mm1 - kk + l;
	t1 = s[k];
	srotg_(&t1, &f, &cs, &sn);
	s[k] = t1;
	if (k == l) {
	    goto L500;
	}
	f = -sn * e[k - 1];
	e[k - 1] = cs * e[k - 1];
L500:
	if (wantv) {
	    srot_(p, &v[k * v_dim1 + 1], &c__1, &v[m * v_dim1 + 1], &c__1, &
		    cs, &sn);
	}
/* L510: */
    }
    goto L610;

/*        SPLIT AT NEGLIGIBLE S(L). */

L520:
    f = e[l - 1];
    e[l - 1] = 0.;
    i__1 = m;
    for (k = l; k <= i__1; ++k) {
	t1 = s[k];
	srotg_(&t1, &f, &cs, &sn);
	s[k] = t1;
	f = -sn * e[k];
	e[k] = cs * e[k];
	if (wantu) {
	    srot_(n, &u[k * u_dim1 + 1], &c__1, &u[(l - 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L530: */
    }
    goto L610;

/*        PERFORM ONE QR STEP. */

L540:

/*           CALCULATE THE SHIFT. */

/* Computing MAX */
    d__6 = (d__1 = s[m], abs(d__1)), d__7 = (d__2 = s[m - 1], abs(d__2)), 
	    d__6 = max(d__6,d__7), d__7 = (d__3 = e[m - 1], abs(d__3)), d__6 =
	     max(d__6,d__7), d__7 = (d__4 = s[l], abs(d__4)), d__6 = max(d__6,
	    d__7), d__7 = (d__5 = e[l], abs(d__5));
    scale = max(d__6,d__7);
    sm = s[m] / scale;
    smm1 = s[m - 1] / scale;
    emm1 = e[m - 1] / scale;
    sl = s[l] / scale;
    el = e[l] / scale;
/* Computing 2nd power */
    d__1 = emm1;
    b = ((smm1 + sm) * (smm1 - sm) + d__1 * d__1) / 2.;
/* Computing 2nd power */
    d__1 = sm * emm1;
    c__ = d__1 * d__1;
    shift = 0.;
    if (b == 0. && c__ == 0.) {
	goto L550;
    }
/* Computing 2nd power */
    d__1 = b;
    shift = sqrt(d__1 * d__1 + c__);
    if (b < 0.) {
	shift = -shift;
    }
    shift = c__ / (b + shift);
L550:
    f = (sl + sm) * (sl - sm) - shift;
    g = sl * el;

/*           CHASE ZEROS. */

    mm1 = m - 1;
    i__1 = mm1;
    for (k = l; k <= i__1; ++k) {
	srotg_(&f, &g, &cs, &sn);
	if (k != l) {
	    e[k - 1] = f;
	}
	f = cs * s[k] + sn * e[k];
	e[k] = cs * e[k] - sn * s[k];
	g = sn * s[k + 1];
	s[k + 1] = cs * s[k + 1];
	if (wantv) {
	    srot_(p, &v[k * v_dim1 + 1], &c__1, &v[(k + 1) * v_dim1 + 1], &
		    c__1, &cs, &sn);
	}
	srotg_(&f, &g, &cs, &sn);
	s[k] = f;
	f = cs * e[k] + sn * s[k + 1];
	s[k + 1] = -sn * e[k] + cs * s[k + 1];
	g = sn * e[k + 1];
	e[k + 1] = cs * e[k + 1];
	if (wantu && k < *n) {
	    srot_(n, &u[k * u_dim1 + 1], &c__1, &u[(k + 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L560: */
    }
    e[m - 1] = f;
    ++iter;
    goto L610;

/*        CONVERGENCE. */

L570:

/*           MAKE THE SINGULAR VALUE  POSITIVE. */

    if (s[l] >= 0.) {
	goto L580;
    }
    s[l] = -s[l];
    if (wantv) {
	sscal_(p, &c_b1033, &v[l * v_dim1 + 1], &c__1);
    }
L580:

/*           ORDER THE SINGULAR VALUE. */

L590:
    if (l == mm) {
	goto L600;
    }
/*           ...EXIT */
    if (s[l] >= s[l + 1]) {
	goto L600;
    }
    t = s[l];
    s[l] = s[l + 1];
    s[l + 1] = t;
    if (wantv && l < *p) {
	sswap_(p, &v[l * v_dim1 + 1], &c__1, &v[(l + 1) * v_dim1 + 1], &c__1);
    }
    if (wantu && l < *n) {
	sswap_(n, &u[l * u_dim1 + 1], &c__1, &u[(l + 1) * u_dim1 + 1], &c__1);
    }
    ++l;
    goto L590;
L600:
    iter = 0;
    --m;
L610:
    goto L360;
L620:
    return 0;
} /* ssvdc_ */

/* Subroutine */ int cgeco_(doublecomplex *a, integer *lda, integer *n, 
	integer *ipvt, doublereal *rcond, doublecomplex *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, l;
    static doublereal s;
    static doublecomplex t;
    static integer kb;
    static doublecomplex ek;
    static doublereal sm;
    static doublecomplex wk;
    static integer kp1;
    static doublecomplex wkm;
    static integer info;
    extern /* Subroutine */ int cgefa_(doublecomplex *, integer *, integer *, 
	    integer *, integer *);
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal anorm;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);


/*     CGECO FACTORS A COMPLEX MATRIX BY GAUSSIAN ELIMINATION */
/*     AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, CGEFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW CGECO BY CGESL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CGECO BY CGESL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CGECO BY CGEDI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW CGECO BY CGEDI. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA, N) */
/*                THE MATRIX TO BE FACTORED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS */
/*                WHICH WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE */
/*                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER */
/*                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR. */

/*        IPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        RCOND   REAL */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       COMPLEX(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK CGEFA */
/*     BLAS CAXPY,CDOTC,CSSCAL,SCASUM */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,CONJG,DBLE */

/*     INTERNAL VARIABLES */



/*     COMPUTE 1-NORM OF A */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --z__;

    /* Function Body */
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = scasum_(n, &a[j * a_dim1 + 1], &c__1);
	anorm = max(d__1,d__2);
/* L10: */
    }

/*     FACTOR */

    cgefa_(&a[a_offset], lda, n, &ipvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E . */
/*     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE CTRANS(U)*W = E */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L20: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) != 0.) {
	    i__3 = k;
	    z__2.r = -z__[i__3].r, z__2.i = -z__[i__3].i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	    d__8 = (d__5 = z__1.r, abs(d__5)) + (d__6 = d_imag(&z__1), abs(
		    d__6));
	    z__4.r = z__1.r / d__8, z__4.i = z__1.i / d__8;
	    z__3.r = d__7 * z__4.r, z__3.i = d__7 * z__4.i;
	    ek.r = z__3.r, ek.i = z__3.i;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = k + k * a_dim1;
	if ((d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(d__2)) <= 
		(d__3 = a[i__3].r, abs(d__3)) + (d__4 = d_imag(&a[k + k * 
		a_dim1]), abs(d__4))) {
	    goto L30;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = k + k * a_dim1;
	s = ((d__1 = a[i__3].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]
		), abs(d__2))) / ((d__3 = z__1.r, abs(d__3)) + (d__4 = d_imag(
		&z__1), abs(d__4)));
	csscal_(n, &s, &z__[1], &c__1);
	z__2.r = s, z__2.i = 0.;
	z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + 
		z__2.i * ek.r;
	ek.r = z__1.r, ek.i = z__1.i;
L30:
	i__2 = k;
	z__1.r = ek.r - z__[i__2].r, z__1.i = ek.i - z__[i__2].i;
	wk.r = z__1.r, wk.i = z__1.i;
	z__2.r = -ek.r, z__2.i = -ek.i;
	i__2 = k;
	z__1.r = z__2.r - z__[i__2].r, z__1.i = z__2.i - z__[i__2].i;
	wkm.r = z__1.r, wkm.i = z__1.i;
	s = (d__1 = wk.r, abs(d__1)) + (d__2 = d_imag(&wk), abs(d__2));
	sm = (d__1 = wkm.r, abs(d__1)) + (d__2 = d_imag(&wkm), abs(d__2));
	i__2 = k + k * a_dim1;
	if ((d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1])
		, abs(d__2)) == 0.) {
	    goto L40;
	}
	d_cnjg(&z__2, &a[k + k * a_dim1]);
	z_div(&z__1, &wk, &z__2);
	wk.r = z__1.r, wk.i = z__1.i;
	d_cnjg(&z__2, &a[k + k * a_dim1]);
	z_div(&z__1, &wkm, &z__2);
	wkm.r = z__1.r, wkm.i = z__1.i;
	goto L50;
L40:
	wk.r = 1., wk.i = 0.;
	wkm.r = 1., wkm.i = 0.;
L50:
	kp1 = k + 1;
	if (kp1 > *n) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    d_cnjg(&z__4, &a[k + j * a_dim1]);
	    z__3.r = wkm.r * z__4.r - wkm.i * z__4.i, z__3.i = wkm.r * z__4.i 
		    + wkm.i * z__4.r;
	    z__2.r = z__[i__3].r + z__3.r, z__2.i = z__[i__3].i + z__3.i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    sm += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(
		    d__2));
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &a[k + j * a_dim1]);
	    z__2.r = wk.r * z__3.r - wk.i * z__3.i, z__2.i = wk.r * z__3.i + 
		    wk.i * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = j;
	    s += (d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&z__[j]), 
		    abs(d__2));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	z__1.r = wkm.r - wk.r, z__1.i = wkm.i - wk.i;
	t.r = z__1.r, t.i = z__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &a[k + j * a_dim1]);
	    z__2.r = t.r * z__3.r - t.i * z__3.i, z__2.i = t.r * z__3.i + t.i 
		    * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
/* L70: */
	}
L80:
L90:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L100: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*     SOLVE CTRANS(L)*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if (k < *n) {
	    i__2 = k;
	    i__3 = k;
	    i__4 = *n - k;
	    cdotc_(&z__2, &i__4, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1], &
		    c__1);
	    z__1.r = z__[i__3].r + z__2.r, z__1.i = z__[i__3].i + z__2.i;
	    z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	}
	i__2 = k;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= 1.) {
	    goto L110;
	}
	i__2 = k;
	s = 1. / ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), 
		abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	i__2 = l;
	t.r = z__[i__2].r, t.i = z__[i__2].i;
	i__2 = l;
	i__3 = k;
	z__[i__2].r = z__[i__3].r, z__[i__2].i = z__[i__3].i;
	i__2 = k;
	z__[i__2].r = t.r, z__[i__2].i = t.i;
/* L120: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE L*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	i__2 = l;
	t.r = z__[i__2].r, t.i = z__[i__2].i;
	i__2 = l;
	i__3 = k;
	z__[i__2].r = z__[i__3].r, z__[i__2].i = z__[i__3].i;
	i__2 = k;
	z__[i__2].r = t.r, z__[i__2].i = t.i;
	if (k < *n) {
	    i__2 = *n - k;
	    caxpy_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1], &
		    c__1);
	}
	i__2 = k;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= 1.) {
	    goto L130;
	}
	i__2 = k;
	s = 1. / ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), 
		abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE  U*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = k + k * a_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= (d__3 = a[i__3].r, abs(d__3)) + (d__4 = d_imag(&a[k 
		+ k * a_dim1]), abs(d__4))) {
	    goto L150;
	}
	i__2 = k + k * a_dim1;
	i__3 = k;
	s = ((d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]
		), abs(d__2))) / ((d__3 = z__[i__3].r, abs(d__3)) + (d__4 = 
		d_imag(&z__[k]), abs(d__4)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	i__2 = k + k * a_dim1;
	if ((d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1])
		, abs(d__2)) != 0.) {
	    i__3 = k;
	    z_div(&z__1, &z__[k], &a[k + k * a_dim1]);
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	}
	i__2 = k + k * a_dim1;
	if ((d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1])
		, abs(d__2)) == 0.) {
	    i__3 = k;
	    z__[i__3].r = 1., z__[i__3].i = 0.;
	}
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
/* L160: */
    }
/*     MAKE ZNORM = 1.0 */
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* cgeco_ */

/* Subroutine */ int cgefa_(doublecomplex *a, integer *lda, integer *n, 
	integer *ipvt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, l;
    static doublecomplex t;
    static integer kp1, nm1;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern integer icamax_(integer *, doublecomplex *, integer *);


/*     CGEFA FACTORS A COMPLEX MATRIX BY GAUSSIAN ELIMINATION. */

/*     CGEFA IS USUALLY CALLED BY CGECO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */
/*     (TIME FOR CGECO) = (1 + 9/N)*(TIME FOR CGEFA) . */

/*     ON ENTRY */

/*        A       COMPLEX(LDA, N) */
/*                THE MATRIX TO BE FACTORED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS */
/*                WHICH WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE */
/*                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER */
/*                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR. */

/*        IPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR */
/*                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES */
/*                     INDICATE THAT CGESL OR CGEDI WILL DIVIDE BY ZERO */
/*                     IF CALLED.  USE  RCOND  IN CGECO FOR A RELIABLE */
/*                     INDICATION OF SINGULARITY. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CSCAL,ICAMAX */
/*     FORTRAN ABS,DIMAG,DBLE */

/*     INTERNAL VARIABLES */



/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;

    /* Function Body */
    *info = 0;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        FIND L = PIVOT INDEX */

	i__2 = *n - k + 1;
	l = icamax_(&i__2, &a[k + k * a_dim1], &c__1) + k - 1;
	ipvt[k] = l;

/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */

	i__2 = l + k * a_dim1;
	if ((d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&a[l + k * a_dim1])
		, abs(d__2)) == 0.) {
	    goto L40;
	}

/*           INTERCHANGE IF NECESSARY */

	if (l == k) {
	    goto L10;
	}
	i__2 = l + k * a_dim1;
	t.r = a[i__2].r, t.i = a[i__2].i;
	i__2 = l + k * a_dim1;
	i__3 = k + k * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = k + k * a_dim1;
	a[i__2].r = t.r, a[i__2].i = t.i;
L10:

/*           COMPUTE MULTIPLIERS */

	z_div(&z__1, &c_b1136, &a[k + k * a_dim1]);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = *n - k;
	cscal_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1);

/*           ROW ELIMINATION WITH COLUMN INDEXING */

	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = l + j * a_dim1;
	    t.r = a[i__3].r, t.i = a[i__3].i;
	    if (l == k) {
		goto L20;
	    }
	    i__3 = l + j * a_dim1;
	    i__4 = k + j * a_dim1;
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
	    i__3 = k + j * a_dim1;
	    a[i__3].r = t.r, a[i__3].i = t.i;
L20:
	    i__3 = *n - k;
	    caxpy_(&i__3, &t, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 + j * 
		    a_dim1], &c__1);
/* L30: */
	}
	goto L50;
L40:
	*info = k;
L50:
/* L60: */
	;
    }
L70:
    ipvt[*n] = *n;
    i__1 = *n + *n * a_dim1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[*n + *n * a_dim1]), 
	    abs(d__2)) == 0.) {
	*info = *n;
    }
    return 0;
} /* cgefa_ */

/* Subroutine */ int cgesl_(doublecomplex *a, integer *lda, integer *n, 
	integer *ipvt, doublecomplex *b, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k, l;
    static doublecomplex t;
    static integer kb, nm1;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CGESL SOLVES THE COMPLEX SYSTEM */
/*     A * X = B  OR  CTRANS(A) * X = B */
/*     USING THE FACTORS COMPUTED BY CGECO OR CGEFA. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA, N) */
/*                THE OUTPUT FROM CGECO OR CGEFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        IPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CGECO OR CGEFA. */

/*        B       COMPLEX(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*        JOB     INTEGER */
/*                = 0         TO SOLVE  A*X = B , */
/*                = NONZERO   TO SOLVE  CTRANS(A)*X = B  WHERE */
/*                            CTRANS(A)  IS THE CONJUGATE TRANSPOSE. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A */
/*        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY */
/*        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER */
/*        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE */
/*        CALLED CORRECTLY AND IF CGECO HAS SET RCOND .GT. 0.0 */
/*        OR CGEFA HAS SET INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL CGECO(A,LDA,N,IPVT,RCOND,Z) */
/*           IF (RCOND IS TOO SMALL) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CGESL(A,LDA,N,IPVT,C(1,J),0) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CDOTC */
/*     FORTRAN CONJG */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --b;

    /* Function Body */
    nm1 = *n - 1;
    if (*job != 0) {
	goto L50;
    }

/*        JOB = 0 , SOLVE  A * X = B */
/*        FIRST SOLVE  L*Y = B */

    if (nm1 < 1) {
	goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	i__2 = l;
	t.r = b[i__2].r, t.i = b[i__2].i;
	if (l == k) {
	    goto L10;
	}
	i__2 = l;
	i__3 = k;
	b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
	i__2 = k;
	b[i__2].r = t.r, b[i__2].i = t.i;
L10:
	i__2 = *n - k;
	caxpy_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &b[k + 1], &c__1);
/* L20: */
    }
L30:

/*        NOW SOLVE  U*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	z_div(&z__1, &b[k], &a[k + k * a_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	i__2 = k;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L40: */
    }
    goto L100;
L50:

/*        JOB = NONZERO, SOLVE  CTRANS(A) * X = B */
/*        FIRST SOLVE  CTRANS(U)*Y = B */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	cdotc_(&z__1, &i__2, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k;
	i__3 = k;
	z__2.r = b[i__3].r - t.r, z__2.i = b[i__3].i - t.i;
	d_cnjg(&z__3, &a[k + k * a_dim1]);
	z_div(&z__1, &z__2, &z__3);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L60: */
    }

/*        NOW SOLVE CTRANS(L)*X = Y */

    if (nm1 < 1) {
	goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
	i__2 = k;
	i__3 = k;
	i__4 = *n - k;
	cdotc_(&z__2, &i__4, &a[k + 1 + k * a_dim1], &c__1, &b[k + 1], &c__1);
	z__1.r = b[i__3].r + z__2.r, z__1.i = b[i__3].i + z__2.i;
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	i__2 = l;
	t.r = b[i__2].r, t.i = b[i__2].i;
	i__2 = l;
	i__3 = k;
	b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
	i__2 = k;
	b[i__2].r = t.r, b[i__2].i = t.i;
L70:
/* L80: */
	;
    }
L90:
L100:
    return 0;
} /* cgesl_ */

/* Subroutine */ int cgedi_(doublecomplex *a, integer *lda, integer *n, 
	integer *ipvt, doublecomplex *det, doublecomplex *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, l;
    static doublecomplex t;
    static integer kb, kp1, nm1;
    static doublereal ten;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), caxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);


/*     CGEDI COMPUTES THE DETERMINANT AND INVERSE OF A MATRIX */
/*     USING THE FACTORS COMPUTED BY CGECO OR CGEFA. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA, N) */
/*                THE OUTPUT FROM CGECO OR CGEFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        IPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CGECO OR CGEFA. */

/*        WORK    COMPLEX(N) */
/*                WORK VECTOR.  CONTENTS DESTROYED. */

/*        JOB     INTEGER */
/*                = 11   BOTH DETERMINANT AND INVERSE. */
/*                = 01   INVERSE ONLY. */
/*                = 10   DETERMINANT ONLY. */

/*     ON RETURN */

/*        A       INVERSE OF ORIGINAL MATRIX IF REQUESTED. */
/*                OTHERWISE UNCHANGED. */

/*        DET     COMPLEX(2) */
/*                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED. */
/*                OTHERWISE NOT REFERENCED. */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. CABS1(DET(1)) .LT. 10.0 */
/*                OR  DET(1) .EQ. 0.0 . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED. */
/*        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY */
/*        AND IF CGECO HAS SET RCOND .GT. 0.0 OR CGEFA HAS SET */
/*        INFO .EQ. 0 . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CSCAL,CSWAP */
/*     FORTRAN ABS,DIMAG,DCMPLX,MOD,DBLE */

/*     INTERNAL VARIABLES */



/*     COMPUTE DETERMINANT */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --det;
    --work;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1].r = 1., det[1].i = 0.;
    det[2].r = 0., det[2].i = 0.;
    ten = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    z__1.r = -det[1].r, z__1.i = -det[1].i;
	    det[1].r = z__1.r, det[1].i = z__1.i;
	}
	i__2 = i__ + i__ * a_dim1;
	z__1.r = a[i__2].r * det[1].r - a[i__2].i * det[1].i, z__1.i = a[i__2]
		.r * det[1].i + a[i__2].i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
/*        ...EXIT */
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 == 0.) {
	    goto L60;
	}
L10:
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 >= 1.) {
	    goto L20;
	}
	z__2.r = ten, z__2.i = 0.;
	z__1.r = z__2.r * det[1].r - z__2.i * det[1].i, z__1.i = z__2.r * det[
		1].i + z__2.i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r - 1., z__1.i = det[2].i - 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L10;
L20:
L30:
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 < ten) {
	    goto L40;
	}
	z__2.r = ten, z__2.i = 0.;
	z_div(&z__1, &det[1], &z__2);
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r + 1., z__1.i = det[2].i + 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*     COMPUTE INVERSE(U) */

    if (*job % 10 == 0) {
	goto L150;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * a_dim1;
	z_div(&z__1, &c_b1092, &a[k + k * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = k + k * a_dim1;
	z__1.r = -a[i__2].r, z__1.i = -a[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	cscal_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = k + j * a_dim1;
	    t.r = a[i__3].r, t.i = a[i__3].i;
	    i__3 = k + j * a_dim1;
	    a[i__3].r = 0., a[i__3].i = 0.;
	    caxpy_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        FORM INVERSE(U)*INVERSE(L) */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L140;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
	kp1 = k + 1;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__ + k * a_dim1;
	    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
	    i__3 = i__ + k * a_dim1;
	    a[i__3].r = 0., a[i__3].i = 0.;
/* L110: */
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    t.r = work[i__3].r, t.i = work[i__3].i;
	    caxpy_(n, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L120: */
	}
	l = ipvt[k];
	if (l != k) {
	    cswap_(n, &a[k * a_dim1 + 1], &c__1, &a[l * a_dim1 + 1], &c__1);
	}
/* L130: */
    }
L140:
L150:
    return 0;
} /* cgedi_ */

/* Subroutine */ int cgbco_(doublecomplex *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublereal *rcond, 
	doublecomplex *z__)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, l, m;
    static doublereal s;
    static doublecomplex t;
    static integer kb, la;
    static doublecomplex ek;
    static integer lm, mm, is, ju;
    static doublereal sm;
    static doublecomplex wk;
    static integer lz, kp1;
    static doublecomplex wkm;
    static integer info;
    extern /* Subroutine */ int cgbfa_(doublecomplex *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal anorm;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);


/*     CGBCO FACTORS A COMPLEX BAND MATRIX BY GAUSSIAN */
/*     ELIMINATION AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, CGBFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW CGBCO BY CGBSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CGBCO BY CGBSL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CGBCO BY CGBDI. */

/*     ON ENTRY */

/*        ABD     COMPLEX(LDA, N) */
/*                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS */
/*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND */
/*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
/*                ML+1 THROUGH 2*ML+MU+1 OF  ABD . */
/*                SEE THE COMMENTS BELOW FOR DETAILS. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */
/*                LDA MUST BE .GE. 2*ML + MU + 1 . */

/*        N       INTEGER */
/*                THE ORDER OF THE ORIGINAL MATRIX. */

/*        ML      INTEGER */
/*                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL. */
/*                0 .LE. ML .LT. N . */

/*        MU      INTEGER */
/*                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */
/*                0 .LE. MU .LT. N . */
/*                MORE EFFICIENT IF  ML .LE. MU . */

/*     ON RETURN */

/*        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
/*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE */
/*                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER */
/*                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR. */

/*        IPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        RCOND   REAL */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       COMPLEX(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     BAND STORAGE */

/*           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT */
/*           WILL SET UP THE INPUT. */

/*                   ML = (BAND WIDTH BELOW THE DIAGONAL) */
/*                   MU = (BAND WIDTH ABOVE THE DIAGONAL) */
/*                   M = ML + MU + 1 */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX0(1, J-MU) */
/*                      I2 = MIN0(N, J+ML) */
/*                      DO 10 I = I1, I2 */
/*                         K = I - J + M */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD . */
/*           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR */
/*           ELEMENTS GENERATED DURING THE TRIANGULARIZATION. */
/*           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 . */
/*           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE */
/*           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED. */

/*     EXAMPLE..  IF THE ORIGINAL MATRIX IS */

/*           11 12 13  0  0  0 */
/*           21 22 23 24  0  0 */
/*            0 32 33 34 35  0 */
/*            0  0 43 44 45 46 */
/*            0  0  0 54 55 56 */
/*            0  0  0  0 65 66 */

/*      THEN  N = 6, ML = 1, MU = 2, LDA .GE. 5  AND ABD SHOULD CONTAIN */

/*            *  *  *  +  +  +  , * = NOT USED */
/*            *  * 13 24 35 46  , + = USED FOR PIVOTING */
/*            * 12 23 34 45 56 */
/*           11 22 33 44 55 66 */
/*           21 32 43 54 65  * */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK CGBFA */
/*     BLAS CAXPY,CDOTC,CSSCAL,SCASUM */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,CONJG,MAX0,MIN0,DBLE */

/*     INTERNAL VARIABLES */



/*     COMPUTE 1-NORM OF A */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --z__;

    /* Function Body */
    anorm = 0.;
    l = *ml + 1;
    is = l + *mu;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = scasum_(&l, &abd[is + j * abd_dim1], &c__1);
	anorm = max(d__1,d__2);
	if (is > *ml + 1) {
	    --is;
	}
	if (j <= *mu) {
	    ++l;
	}
	if (j >= *n - *ml) {
	    --l;
	}
/* L10: */
    }

/*     FACTOR */

    cgbfa_(&abd[abd_offset], lda, n, ml, mu, &ipvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  CTRANS(A)*Y = E . */
/*     CTRANS(A)  IS THE CONJUGATE TRANSPOSE OF A . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(U)*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE CTRANS(U)*W = E */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L20: */
    }
    m = *ml + *mu + 1;
    ju = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) != 0.) {
	    i__3 = k;
	    z__2.r = -z__[i__3].r, z__2.i = -z__[i__3].i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	    d__8 = (d__5 = z__1.r, abs(d__5)) + (d__6 = d_imag(&z__1), abs(
		    d__6));
	    z__4.r = z__1.r / d__8, z__4.i = z__1.i / d__8;
	    z__3.r = d__7 * z__4.r, z__3.i = d__7 * z__4.i;
	    ek.r = z__3.r, ek.i = z__3.i;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = m + k * abd_dim1;
	if ((d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(d__2)) <= 
		(d__3 = abd[i__3].r, abs(d__3)) + (d__4 = d_imag(&abd[m + k * 
		abd_dim1]), abs(d__4))) {
	    goto L30;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = m + k * abd_dim1;
	s = ((d__1 = abd[i__3].r, abs(d__1)) + (d__2 = d_imag(&abd[m + k * 
		abd_dim1]), abs(d__2))) / ((d__3 = z__1.r, abs(d__3)) + (d__4 
		= d_imag(&z__1), abs(d__4)));
	csscal_(n, &s, &z__[1], &c__1);
	z__2.r = s, z__2.i = 0.;
	z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + 
		z__2.i * ek.r;
	ek.r = z__1.r, ek.i = z__1.i;
L30:
	i__2 = k;
	z__1.r = ek.r - z__[i__2].r, z__1.i = ek.i - z__[i__2].i;
	wk.r = z__1.r, wk.i = z__1.i;
	z__2.r = -ek.r, z__2.i = -ek.i;
	i__2 = k;
	z__1.r = z__2.r - z__[i__2].r, z__1.i = z__2.i - z__[i__2].i;
	wkm.r = z__1.r, wkm.i = z__1.i;
	s = (d__1 = wk.r, abs(d__1)) + (d__2 = d_imag(&wk), abs(d__2));
	sm = (d__1 = wkm.r, abs(d__1)) + (d__2 = d_imag(&wkm), abs(d__2));
	i__2 = m + k * abd_dim1;
	if ((d__1 = abd[i__2].r, abs(d__1)) + (d__2 = d_imag(&abd[m + k * 
		abd_dim1]), abs(d__2)) == 0.) {
	    goto L40;
	}
	d_cnjg(&z__2, &abd[m + k * abd_dim1]);
	z_div(&z__1, &wk, &z__2);
	wk.r = z__1.r, wk.i = z__1.i;
	d_cnjg(&z__2, &abd[m + k * abd_dim1]);
	z_div(&z__1, &wkm, &z__2);
	wkm.r = z__1.r, wkm.i = z__1.i;
	goto L50;
L40:
	wk.r = 1., wk.i = 0.;
	wkm.r = 1., wkm.i = 0.;
L50:
	kp1 = k + 1;
/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = max(i__3,i__4);
	ju = min(i__2,*n);
	mm = m;
	if (kp1 > ju) {
	    goto L90;
	}
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --mm;
	    i__3 = j;
	    d_cnjg(&z__4, &abd[mm + j * abd_dim1]);
	    z__3.r = wkm.r * z__4.r - wkm.i * z__4.i, z__3.i = wkm.r * z__4.i 
		    + wkm.i * z__4.r;
	    z__2.r = z__[i__3].r + z__3.r, z__2.i = z__[i__3].i + z__3.i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    sm += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(
		    d__2));
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &abd[mm + j * abd_dim1]);
	    z__2.r = wk.r * z__3.r - wk.i * z__3.i, z__2.i = wk.r * z__3.i + 
		    wk.i * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = j;
	    s += (d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&z__[j]), 
		    abs(d__2));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	z__1.r = wkm.r - wk.r, z__1.i = wkm.i - wk.i;
	t.r = z__1.r, t.i = z__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	mm = m;
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --mm;
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &abd[mm + j * abd_dim1]);
	    z__2.r = t.r * z__3.r - t.i * z__3.i, z__2.i = t.r * z__3.i + t.i 
		    * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
/* L70: */
	}
L80:
L90:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L100: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*     SOLVE CTRANS(L)*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	if (k < *n) {
	    i__2 = k;
	    i__3 = k;
	    cdotc_(&z__2, &lm, &abd[m + 1 + k * abd_dim1], &c__1, &z__[k + 1],
		     &c__1);
	    z__1.r = z__[i__3].r + z__2.r, z__1.i = z__[i__3].i + z__2.i;
	    z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	}
	i__2 = k;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= 1.) {
	    goto L110;
	}
	i__2 = k;
	s = 1. / ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), 
		abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	i__2 = l;
	t.r = z__[i__2].r, t.i = z__[i__2].i;
	i__2 = l;
	i__3 = k;
	z__[i__2].r = z__[i__3].r, z__[i__2].i = z__[i__3].i;
	i__2 = k;
	z__[i__2].r = t.r, z__[i__2].i = t.i;
/* L120: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE L*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	i__2 = l;
	t.r = z__[i__2].r, t.i = z__[i__2].i;
	i__2 = l;
	i__3 = k;
	z__[i__2].r = z__[i__3].r, z__[i__2].i = z__[i__3].i;
	i__2 = k;
	z__[i__2].r = t.r, z__[i__2].i = t.i;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	if (k < *n) {
	    caxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &z__[k + 1], &
		    c__1);
	}
	i__2 = k;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= 1.) {
	    goto L130;
	}
	i__2 = k;
	s = 1. / ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), 
		abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE  U*Z = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = m + k * abd_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= (d__3 = abd[i__3].r, abs(d__3)) + (d__4 = d_imag(&
		abd[m + k * abd_dim1]), abs(d__4))) {
	    goto L150;
	}
	i__2 = m + k * abd_dim1;
	i__3 = k;
	s = ((d__1 = abd[i__2].r, abs(d__1)) + (d__2 = d_imag(&abd[m + k * 
		abd_dim1]), abs(d__2))) / ((d__3 = z__[i__3].r, abs(d__3)) + (
		d__4 = d_imag(&z__[k]), abs(d__4)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	i__2 = m + k * abd_dim1;
	if ((d__1 = abd[i__2].r, abs(d__1)) + (d__2 = d_imag(&abd[m + k * 
		abd_dim1]), abs(d__2)) != 0.) {
	    i__3 = k;
	    z_div(&z__1, &z__[k], &abd[m + k * abd_dim1]);
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	}
	i__2 = m + k * abd_dim1;
	if ((d__1 = abd[i__2].r, abs(d__1)) + (d__2 = d_imag(&abd[m + k * 
		abd_dim1]), abs(d__2)) == 0.) {
	    i__3 = k;
	    z__[i__3].r = 1., z__[i__3].i = 0.;
	}
	lm = min(k,m) - 1;
	la = m - lm;
	lz = k - lm;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	caxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lz], &c__1);
/* L160: */
    }
/*     MAKE ZNORM = 1.0 */
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* cgbco_ */

/* Subroutine */ int cgbfa_(doublecomplex *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublecomplex t;
    static integer i0, j0, j1, lm, mm, ju, jz, kp1, nm1;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern integer icamax_(integer *, doublecomplex *, integer *);


/*     CGBFA FACTORS A COMPLEX BAND MATRIX BY ELIMINATION. */

/*     CGBFA IS USUALLY CALLED BY CGBCO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */

/*     ON ENTRY */

/*        ABD     COMPLEX(LDA, N) */
/*                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS */
/*                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND */
/*                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS */
/*                ML+1 THROUGH 2*ML+MU+1 OF  ABD . */
/*                SEE THE COMMENTS BELOW FOR DETAILS. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */
/*                LDA MUST BE .GE. 2*ML + MU + 1 . */

/*        N       INTEGER */
/*                THE ORDER OF THE ORIGINAL MATRIX. */

/*        ML      INTEGER */
/*                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL. */
/*                0 .LE. ML .LT. N . */

/*        MU      INTEGER */
/*                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */
/*                0 .LE. MU .LT. N . */
/*                MORE EFFICIENT IF  ML .LE. MU . */
/*     ON RETURN */

/*        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND */
/*                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = L*U  WHERE */
/*                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER */
/*                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR. */

/*        IPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR */
/*                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES */
/*                     INDICATE THAT CGBSL WILL DIVIDE BY ZERO IF */
/*                     CALLED.  USE  RCOND  IN CGBCO FOR A RELIABLE */
/*                     INDICATION OF SINGULARITY. */

/*     BAND STORAGE */

/*           IF  A  IS A BAND MATRIX, THE FOLLOWING PROGRAM SEGMENT */
/*           WILL SET UP THE INPUT. */

/*                   ML = (BAND WIDTH BELOW THE DIAGONAL) */
/*                   MU = (BAND WIDTH ABOVE THE DIAGONAL) */
/*                   M = ML + MU + 1 */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX0(1, J-MU) */
/*                      I2 = MIN0(N, J+ML) */
/*                      DO 10 I = I1, I2 */
/*                         K = I - J + M */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD . */
/*           IN ADDITION, THE FIRST  ML  ROWS IN  ABD  ARE USED FOR */
/*           ELEMENTS GENERATED DURING THE TRIANGULARIZATION. */
/*           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 . */
/*           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE */
/*           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CSCAL,ICAMAX */
/*     FORTRAN ABS,DIMAG,MAX0,MIN0,DBLE */

/*     INTERNAL VARIABLES */



    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;

    /* Function Body */
    m = *ml + *mu + 1;
    *info = 0;

/*     ZERO INITIAL FILL-IN COLUMNS */

    j0 = *mu + 2;
    j1 = min(*n,m) - 1;
    if (j1 < j0) {
	goto L30;
    }
    i__1 = j1;
    for (jz = j0; jz <= i__1; ++jz) {
	i0 = m + 1 - jz;
	i__2 = *ml;
	for (i__ = i0; i__ <= i__2; ++i__) {
	    i__3 = i__ + jz * abd_dim1;
	    abd[i__3].r = 0., abd[i__3].i = 0.;
/* L10: */
	}
/* L20: */
    }
L30:
    jz = j1;
    ju = 0;

/*     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L130;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        ZERO NEXT FILL-IN COLUMN */

	++jz;
	if (jz > *n) {
	    goto L50;
	}
	if (*ml < 1) {
	    goto L50;
	}
	i__2 = *ml;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + jz * abd_dim1;
	    abd[i__3].r = 0., abd[i__3].i = 0.;
/* L40: */
	}
L50:

/*        FIND L = PIVOT INDEX */

/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	i__2 = lm + 1;
	l = icamax_(&i__2, &abd[m + k * abd_dim1], &c__1) + m - 1;
	ipvt[k] = l + k - m;

/*        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED */

	i__2 = l + k * abd_dim1;
	if ((d__1 = abd[i__2].r, abs(d__1)) + (d__2 = d_imag(&abd[l + k * 
		abd_dim1]), abs(d__2)) == 0.) {
	    goto L100;
	}

/*           INTERCHANGE IF NECESSARY */

	if (l == m) {
	    goto L60;
	}
	i__2 = l + k * abd_dim1;
	t.r = abd[i__2].r, t.i = abd[i__2].i;
	i__2 = l + k * abd_dim1;
	i__3 = m + k * abd_dim1;
	abd[i__2].r = abd[i__3].r, abd[i__2].i = abd[i__3].i;
	i__2 = m + k * abd_dim1;
	abd[i__2].r = t.r, abd[i__2].i = t.i;
L60:

/*           COMPUTE MULTIPLIERS */

	z_div(&z__1, &c_b1136, &abd[m + k * abd_dim1]);
	t.r = z__1.r, t.i = z__1.i;
	cscal_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1);

/*           ROW ELIMINATION WITH COLUMN INDEXING */

/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = max(i__3,i__4);
	ju = min(i__2,*n);
	mm = m;
	if (ju < kp1) {
	    goto L90;
	}
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --l;
	    --mm;
	    i__3 = l + j * abd_dim1;
	    t.r = abd[i__3].r, t.i = abd[i__3].i;
	    if (l == mm) {
		goto L70;
	    }
	    i__3 = l + j * abd_dim1;
	    i__4 = mm + j * abd_dim1;
	    abd[i__3].r = abd[i__4].r, abd[i__3].i = abd[i__4].i;
	    i__3 = mm + j * abd_dim1;
	    abd[i__3].r = t.r, abd[i__3].i = t.i;
L70:
	    caxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &abd[mm + 1 + 
		    j * abd_dim1], &c__1);
/* L80: */
	}
L90:
	goto L110;
L100:
	*info = k;
L110:
/* L120: */
	;
    }
L130:
    ipvt[*n] = *n;
    i__1 = m + *n * abd_dim1;
    if ((d__1 = abd[i__1].r, abs(d__1)) + (d__2 = d_imag(&abd[m + *n * 
	    abd_dim1]), abs(d__2)) == 0.) {
	*info = *n;
    }
    return 0;
} /* cgbfa_ */

/* Subroutine */ int cgbsl_(doublecomplex *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublecomplex *b, integer *
	job)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k, l, m;
    static doublecomplex t;
    static integer kb, la, lb, lm, nm1;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CGBSL SOLVES THE COMPLEX BAND SYSTEM */
/*     A * X = B  OR  CTRANS(A) * X = B */
/*     USING THE FACTORS COMPUTED BY CGBCO OR CGBFA. */

/*     ON ENTRY */

/*        ABD     COMPLEX(LDA, N) */
/*                THE OUTPUT FROM CGBCO OR CGBFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */

/*        N       INTEGER */
/*                THE ORDER OF THE ORIGINAL MATRIX. */

/*        ML      INTEGER */
/*                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL. */

/*        MU      INTEGER */
/*                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */

/*        IPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CGBCO OR CGBFA. */

/*        B       COMPLEX(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*        JOB     INTEGER */
/*                = 0         TO SOLVE  A*X = B , */
/*                = NONZERO   TO SOLVE  CTRANS(A)*X = B , WHERE */
/*                            CTRANS(A)  IS THE CONJUGATE TRANSPOSE. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A */
/*        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY */
/*        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER */
/*        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE */
/*        CALLED CORRECTLY AND IF CGBCO HAS SET RCOND .GT. 0.0 */
/*        OR CGBFA HAS SET INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL CGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z) */
/*           IF (RCOND IS TOO SMALL) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CDOTC */
/*     FORTRAN CONJG,MIN0 */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --b;

    /* Function Body */
    m = *mu + *ml + 1;
    nm1 = *n - 1;
    if (*job != 0) {
	goto L50;
    }

/*        JOB = 0 , SOLVE  A * X = B */
/*        FIRST SOLVE L*Y = B */

    if (*ml == 0) {
	goto L30;
    }
    if (nm1 < 1) {
	goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	l = ipvt[k];
	i__2 = l;
	t.r = b[i__2].r, t.i = b[i__2].i;
	if (l == k) {
	    goto L10;
	}
	i__2 = l;
	i__3 = k;
	b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
	i__2 = k;
	b[i__2].r = t.r, b[i__2].i = t.i;
L10:
	caxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &c__1);
/* L20: */
    }
L30:

/*        NOW SOLVE  U*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	z_div(&z__1, &b[k], &abd[m + k * abd_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	lm = min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	i__2 = k;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	caxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/* L40: */
    }
    goto L100;
L50:

/*        JOB = NONZERO, SOLVE  CTRANS(A) * X = B */
/*        FIRST SOLVE  CTRANS(U)*Y = B */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	lm = min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	cdotc_(&z__1, &lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k;
	i__3 = k;
	z__2.r = b[i__3].r - t.r, z__2.i = b[i__3].i - t.i;
	d_cnjg(&z__3, &abd[m + k * abd_dim1]);
	z_div(&z__1, &z__2, &z__3);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L60: */
    }

/*        NOW SOLVE CTRANS(L)*X = Y */

    if (*ml == 0) {
	goto L90;
    }
    if (nm1 < 1) {
	goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	i__2 = k;
	i__3 = k;
	cdotc_(&z__2, &lm, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &
		c__1);
	z__1.r = b[i__3].r + z__2.r, z__1.i = b[i__3].i + z__2.i;
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	i__2 = l;
	t.r = b[i__2].r, t.i = b[i__2].i;
	i__2 = l;
	i__3 = k;
	b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
	i__2 = k;
	b[i__2].r = t.r, b[i__2].i = t.i;
L70:
/* L80: */
	;
    }
L90:
L100:
    return 0;
} /* cgbsl_ */

/* Subroutine */ int cgbdi_(doublecomplex *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublecomplex *det)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, m;
    static doublereal ten;


/*     CGBDI COMPUTES THE DETERMINANT OF A BAND MATRIX */
/*     USING THE FACTORS COMPUTED BY CGBCO OR CGBFA. */
/*     IF THE INVERSE IS NEEDED, USE CGBSL  N  TIMES. */

/*     ON ENTRY */

/*        ABD     COMPLEX(LDA, N) */
/*                THE OUTPUT FROM CGBCO OR CGBFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */

/*        N       INTEGER */
/*                THE ORDER OF THE ORIGINAL MATRIX. */

/*        ML      INTEGER */
/*                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL. */

/*        MU      INTEGER */
/*                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */

/*        IPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CGBCO OR CGBFA. */

/*     ON RETURN */

/*        DET     COMPLEX(2) */
/*                DETERMINANT OF ORIGINAL MATRIX. */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. CABS1(DET(1)) .LT. 10.0 */
/*                OR  DET(1) = 0.0 . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     FORTRAN ABS,DIMAG,DCMPLX,DBLE */

/*     INTERNAL VARIABLES */



    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --det;

    /* Function Body */
    m = *ml + *mu + 1;
    det[1].r = 1., det[1].i = 0.;
    det[2].r = 0., det[2].i = 0.;
    ten = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    z__1.r = -det[1].r, z__1.i = -det[1].i;
	    det[1].r = z__1.r, det[1].i = z__1.i;
	}
	i__2 = m + i__ * abd_dim1;
	z__1.r = abd[i__2].r * det[1].r - abd[i__2].i * det[1].i, z__1.i = 
		abd[i__2].r * det[1].i + abd[i__2].i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
/*     ...EXIT */
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 == 0.) {
	    goto L60;
	}
L10:
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 >= 1.) {
	    goto L20;
	}
	z__2.r = ten, z__2.i = 0.;
	z__1.r = z__2.r * det[1].r - z__2.i * det[1].i, z__1.i = z__2.r * det[
		1].i + z__2.i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r - 1., z__1.i = det[2].i - 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L10;
L20:
L30:
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 < ten) {
	    goto L40;
	}
	z__2.r = ten, z__2.i = 0.;
	z_div(&z__1, &det[1], &z__2);
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r + 1., z__1.i = det[2].i + 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
    return 0;
} /* cgbdi_ */

/* Subroutine */ int cpoco_(doublecomplex *a, integer *lda, integer *n, 
	doublereal *rcond, doublecomplex *z__, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static doublecomplex t;
    static integer kb;
    static doublecomplex ek;
    static doublereal sm;
    static doublecomplex wk;
    static integer jm1, kp1;
    static doublecomplex wkm;
    extern /* Subroutine */ int cpofa_(doublecomplex *, integer *, integer *, 
	    integer *);
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal anorm;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);


/*     CPOCO FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX */
/*     AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, CPOFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW CPOCO BY CPOSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CPOCO BY CPOSL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CPOCO BY CPODI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW CPOCO BY CPODI. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA, N) */
/*                THE HERMITIAN MATRIX TO BE FACTORED.  ONLY THE */
/*                DIAGONAL AND UPPER TRIANGLE ARE USED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = */
/*                CTRANS(R)*R WHERE  CTRANS(R)  IS THE CONJUGATE */
/*                TRANSPOSE.  THE STRICT LOWER TRIANGLE IS UNALTERED. */
/*                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE. */

/*        RCOND   REAL */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED. */

/*        Z       COMPLEX(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */
/*                IF  INFO .NE. 0 , Z  IS UNCHANGED. */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR */
/*                     OF ORDER  K  IS NOT POSITIVE DEFINITE. */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK CPOFA */
/*     BLAS CAXPY,CDOTC,CSSCAL,SCASUM */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,CONJG,DBLE */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A USING ONLY UPPER HALF */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	d__1 = scasum_(&j, &a[j * a_dim1 + 1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__ + j * a_dim1;
	    d__3 = z__[i__4].r + ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[i__ + j * a_dim1]), abs(d__2)));
	    z__1.r = d__3, z__1.i = 0.;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = j;
	d__1 = anorm, d__2 = z__[i__2].r;
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    cpofa_(&a[a_offset], lda, n, info);
    if (*info != 0) {
	goto L180;
    }

/*        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*        GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(R)*W = E . */
/*        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*        SOLVE CTRANS(R)*W = E */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L50: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) != 0.) {
	    i__3 = k;
	    z__2.r = -z__[i__3].r, z__2.i = -z__[i__3].i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	    d__8 = (d__5 = z__1.r, abs(d__5)) + (d__6 = d_imag(&z__1), abs(
		    d__6));
	    z__4.r = z__1.r / d__8, z__4.i = z__1.i / d__8;
	    z__3.r = d__7 * z__4.r, z__3.i = d__7 * z__4.i;
	    ek.r = z__3.r, ek.i = z__3.i;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = k + k * a_dim1;
	if ((d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(d__2)) <= 
		a[i__3].r) {
	    goto L60;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = k + k * a_dim1;
	s = a[i__3].r / ((d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
		abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
	z__2.r = s, z__2.i = 0.;
	z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + 
		z__2.i * ek.r;
	ek.r = z__1.r, ek.i = z__1.i;
L60:
	i__2 = k;
	z__1.r = ek.r - z__[i__2].r, z__1.i = ek.i - z__[i__2].i;
	wk.r = z__1.r, wk.i = z__1.i;
	z__2.r = -ek.r, z__2.i = -ek.i;
	i__2 = k;
	z__1.r = z__2.r - z__[i__2].r, z__1.i = z__2.i - z__[i__2].i;
	wkm.r = z__1.r, wkm.i = z__1.i;
	s = (d__1 = wk.r, abs(d__1)) + (d__2 = d_imag(&wk), abs(d__2));
	sm = (d__1 = wkm.r, abs(d__1)) + (d__2 = d_imag(&wkm), abs(d__2));
	z_div(&z__1, &wk, &a[k + k * a_dim1]);
	wk.r = z__1.r, wk.i = z__1.i;
	z_div(&z__1, &wkm, &a[k + k * a_dim1]);
	wkm.r = z__1.r, wkm.i = z__1.i;
	kp1 = k + 1;
	if (kp1 > *n) {
	    goto L100;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    d_cnjg(&z__4, &a[k + j * a_dim1]);
	    z__3.r = wkm.r * z__4.r - wkm.i * z__4.i, z__3.i = wkm.r * z__4.i 
		    + wkm.i * z__4.r;
	    z__2.r = z__[i__3].r + z__3.r, z__2.i = z__[i__3].i + z__3.i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    sm += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(
		    d__2));
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &a[k + j * a_dim1]);
	    z__2.r = wk.r * z__3.r - wk.i * z__3.i, z__2.i = wk.r * z__3.i + 
		    wk.i * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = j;
	    s += (d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&z__[j]), 
		    abs(d__2));
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	z__1.r = wkm.r - wk.r, z__1.i = wkm.i - wk.i;
	t.r = z__1.r, t.i = z__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &a[k + j * a_dim1]);
	    z__2.r = t.r * z__3.r - t.i * z__3.i, z__2.i = t.r * z__3.i + t.i 
		    * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
/* L80: */
	}
L90:
L100:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L110: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*        SOLVE R*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = k + k * a_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= a[i__3].r) {
	    goto L120;
	}
	i__2 = k + k * a_dim1;
	i__3 = k;
	s = a[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		z__[k]), abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
L120:
	i__2 = k;
	z_div(&z__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
/* L130: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*        SOLVE CTRANS(R)*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	i__3 = k;
	i__4 = k - 1;
	cdotc_(&z__2, &i__4, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
	z__1.r = z__[i__3].r - z__2.r, z__1.i = z__[i__3].i - z__2.i;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	i__2 = k;
	i__3 = k + k * a_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= a[i__3].r) {
	    goto L140;
	}
	i__2 = k + k * a_dim1;
	i__3 = k;
	s = a[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		z__[k]), abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	i__2 = k;
	z_div(&z__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
/* L150: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        SOLVE R*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = k + k * a_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= a[i__3].r) {
	    goto L160;
	}
	i__2 = k + k * a_dim1;
	i__3 = k;
	s = a[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		z__[k]), abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	i__2 = k;
	z_div(&z__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
/* L170: */
    }
/*        MAKE ZNORM = 1.0 */
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
L180:
    return 0;
} /* cpoco_ */

/* Subroutine */ int cpofa_(doublecomplex *a, integer *lda, integer *n, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal s;
    static doublecomplex t;
    static integer jm1;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CPOFA FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX. */

/*     CPOFA IS USUALLY CALLED BY CPOCO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */
/*     (TIME FOR CPOCO) = (1 + 18/N)*(TIME FOR CPOFA) . */

/*     ON ENTRY */

/*        A       COMPLEX(LDA, N) */
/*                THE HERMITIAN MATRIX TO BE FACTORED.  ONLY THE */
/*                DIAGONAL AND UPPER TRIANGLE ARE USED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = */
/*                CTRANS(R)*R WHERE  CTRANS(R)  IS THE CONJUGATE */
/*                TRANSPOSE.  THE STRICT LOWER TRIANGLE IS UNALTERED. */
/*                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE. */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR */
/*                     OF ORDER  K  IS NOT POSITIVE DEFINITE. */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CDOTC */
/*     FORTRAN DIMAG,DCMPLX,CONJG,DBLE,SQRT */

/*     INTERNAL VARIABLES */

/*     BEGIN BLOCK WITH ...EXITS TO 40 */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = k + j * a_dim1;
	    i__4 = k - 1;
	    cdotc_(&z__2, &i__4, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1]
		    , &c__1);
	    z__1.r = a[i__3].r - z__2.r, z__1.i = a[i__3].i - z__2.i;
	    t.r = z__1.r, t.i = z__1.i;
	    z_div(&z__1, &t, &a[k + k * a_dim1]);
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = k + j * a_dim1;
	    a[i__3].r = t.r, a[i__3].i = t.i;
	    d_cnjg(&z__2, &t);
	    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i 
		    * z__2.r;
	    s += z__1.r;
/* L10: */
	}
L20:
	i__2 = j + j * a_dim1;
	s = a[i__2].r - s;
/*     ......EXIT */
	if (s <= 0. || d_imag(&a[j + j * a_dim1]) != 0.) {
	    goto L40;
	}
	i__2 = j + j * a_dim1;
	d__1 = sqrt(s);
	z__1.r = d__1, z__1.i = 0.;
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* cpofa_ */

/* Subroutine */ int cposl_(doublecomplex *a, integer *lda, integer *n, 
	doublecomplex *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex t;
    static integer kb;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CPOSL SOLVES THE COMPLEX HERMITIAN POSITIVE DEFINITE SYSTEM */
/*     A * X = B */
/*     USING THE FACTORS COMPUTED BY CPOCO OR CPOFA. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA, N) */
/*                THE OUTPUT FROM CPOCO OR CPOFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        B       COMPLEX(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES */
/*        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE */
/*        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED */
/*        CORRECTLY AND  INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL CPOCO(A,LDA,N,RCOND,Z,INFO) */
/*           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CPOSL(A,LDA,N,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CDOTC */

/*     INTERNAL VARIABLES */


/*     SOLVE CTRANS(R)*Y = B */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	cdotc_(&z__1, &i__2, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k;
	i__3 = k;
	z__2.r = b[i__3].r - t.r, z__2.i = b[i__3].i - t.i;
	z_div(&z__1, &z__2, &a[k + k * a_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L10: */
    }

/*     SOLVE R*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	z_div(&z__1, &b[k], &a[k + k * a_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	i__2 = k;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L20: */
    }
    return 0;
} /* cposl_ */

/* Subroutine */ int cpodi_(doublecomplex *a, integer *lda, integer *n, 
	doublereal *det, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static doublecomplex t;
    static integer jm1, kp1;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CPODI COMPUTES THE DETERMINANT AND INVERSE OF A CERTAIN */
/*     COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX (SEE BELOW) */
/*     USING THE FACTORS COMPUTED BY CPOCO, CPOFA OR CQRDC. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA, N) */
/*                THE OUTPUT  A  FROM CPOCO OR CPOFA */
/*                OR THE OUTPUT  X  FROM CQRDC. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        JOB     INTEGER */
/*                = 11   BOTH DETERMINANT AND INVERSE. */
/*                = 01   INVERSE ONLY. */
/*                = 10   DETERMINANT ONLY. */

/*     ON RETURN */

/*        A       IF CPOCO OR CPOFA WAS USED TO FACTOR  A  THEN */
/*                CPODI PRODUCES THE UPPER HALF OF INVERSE(A) . */
/*                IF CQRDC WAS USED TO DECOMPOSE  X  THEN */
/*                CPODI PRODUCES THE UPPER HALF OF INVERSE(CTRANS(X)*X) */
/*                WHERE CTRANS(X) IS THE CONJUGATE TRANSPOSE. */
/*                ELEMENTS OF  A  BELOW THE DIAGONAL ARE UNCHANGED. */
/*                IF THE UNITS DIGIT OF JOB IS ZERO,  A  IS UNCHANGED. */

/*        DET     REAL(2) */
/*                DETERMINANT OF  A  OR OF  CTRANS(X)*X  IF REQUESTED. */
/*                OTHERWISE NOT REFERENCED. */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. DET(1) .LT. 10.0 */
/*                OR  DET(1) .EQ. 0.0 . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED. */
/*        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY */
/*        AND IF CPOCO OR CPOFA HAS SET INFO .EQ. 0 . */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CSCAL */
/*     FORTRAN CONJG,MOD,DBLE */

/*     INTERNAL VARIABLES */


/*     COMPUTE DETERMINANT */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --det;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    s = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + i__ * a_dim1;
/* Computing 2nd power */
	d__1 = a[i__2].r;
	det[1] = d__1 * d__1 * det[1];
/*        ...EXIT */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (det[1] >= 1.) {
	    goto L20;
	}
	det[1] = s * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (det[1] < s) {
	    goto L40;
	}
	det[1] /= s;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*     COMPUTE INVERSE(R) */

    if (*job % 10 == 0) {
	goto L140;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * a_dim1;
	z_div(&z__1, &c_b1092, &a[k + k * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = k + k * a_dim1;
	z__1.r = -a[i__2].r, z__1.i = -a[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	cscal_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = k + j * a_dim1;
	    t.r = a[i__3].r, t.i = a[i__3].i;
	    i__3 = k + j * a_dim1;
	    a[i__3].r = 0., a[i__3].i = 0.;
	    caxpy_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        FORM  INVERSE(R) * CTRANS(INVERSE(R)) */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L120;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    d_cnjg(&z__1, &a[k + j * a_dim1]);
	    t.r = z__1.r, t.i = z__1.i;
	    caxpy_(&k, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L110: */
	}
L120:
	d_cnjg(&z__1, &a[j + j * a_dim1]);
	t.r = z__1.r, t.i = z__1.i;
	cscal_(&j, &t, &a[j * a_dim1 + 1], &c__1);
/* L130: */
    }
L140:
    return 0;
} /* cpodi_ */

/* Subroutine */ int cppco_(doublecomplex *ap, integer *n, doublereal *rcond, 
	doublecomplex *z__, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static doublecomplex t;
    static integer j1, kb;
    static doublecomplex ek;
    static integer ij, kj, kk;
    static doublereal sm;
    static doublecomplex wk;
    static integer jm1, kp1;
    static doublecomplex wkm;
    extern /* Subroutine */ int cppfa_(doublecomplex *, integer *, integer *);
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal anorm;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);


/*     CPPCO FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX */
/*     STORED IN PACKED FORM */
/*     AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, CPPFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW CPPCO BY CPPSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CPPCO BY CPPSL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CPPCO BY CPPDI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW CPPCO BY CPPDI. */

/*     ON ENTRY */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                THE PACKED FORM OF A HERMITIAN MATRIX  A .  THE */
/*                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY */
/*                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 . */
/*                SEE COMMENTS BELOW FOR DETAILS. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        AP      AN UPPER TRIANGULAR MATRIX  R , STORED IN PACKED */
/*                FORM, SO THAT  A = CTRANS(R)*R . */
/*                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE. */

/*        RCOND   REAL */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED. */

/*        Z       COMPLEX(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */
/*                IF  INFO .NE. 0 , Z  IS UNCHANGED. */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR */
/*                     OF ORDER  K  IS NOT POSITIVE DEFINITE. */

/*     PACKED STORAGE */

/*          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER */
/*          TRIANGLE OF A HERMITIAN MATRIX. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K) = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK CPPFA */
/*     BLAS CAXPY,CDOTC,CSSCAL,SCASUM */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,CONJG,DBLE */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A */

    /* Parameter adjustments */
    --z__;
    --ap;

    /* Function Body */
    j1 = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	d__1 = scasum_(&j, &ap[j1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	ij = j1;
	j1 += j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = ij;
	    d__3 = z__[i__4].r + ((d__1 = ap[i__5].r, abs(d__1)) + (d__2 = 
		    d_imag(&ap[ij]), abs(d__2)));
	    z__1.r = d__3, z__1.i = 0.;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    ++ij;
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = j;
	d__1 = anorm, d__2 = z__[i__2].r;
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    cppfa_(&ap[1], n, info);
    if (*info != 0) {
	goto L180;
    }

/*        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*        GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(R)*W = E . */
/*        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*        SOLVE CTRANS(R)*W = E */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L50: */
    }
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk += k;
	i__2 = k;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) != 0.) {
	    i__3 = k;
	    z__2.r = -z__[i__3].r, z__2.i = -z__[i__3].i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	    d__8 = (d__5 = z__1.r, abs(d__5)) + (d__6 = d_imag(&z__1), abs(
		    d__6));
	    z__4.r = z__1.r / d__8, z__4.i = z__1.i / d__8;
	    z__3.r = d__7 * z__4.r, z__3.i = d__7 * z__4.i;
	    ek.r = z__3.r, ek.i = z__3.i;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = kk;
	if ((d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(d__2)) <= 
		ap[i__3].r) {
	    goto L60;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = kk;
	s = ap[i__3].r / ((d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), 
		abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
	z__2.r = s, z__2.i = 0.;
	z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + 
		z__2.i * ek.r;
	ek.r = z__1.r, ek.i = z__1.i;
L60:
	i__2 = k;
	z__1.r = ek.r - z__[i__2].r, z__1.i = ek.i - z__[i__2].i;
	wk.r = z__1.r, wk.i = z__1.i;
	z__2.r = -ek.r, z__2.i = -ek.i;
	i__2 = k;
	z__1.r = z__2.r - z__[i__2].r, z__1.i = z__2.i - z__[i__2].i;
	wkm.r = z__1.r, wkm.i = z__1.i;
	s = (d__1 = wk.r, abs(d__1)) + (d__2 = d_imag(&wk), abs(d__2));
	sm = (d__1 = wkm.r, abs(d__1)) + (d__2 = d_imag(&wkm), abs(d__2));
	z_div(&z__1, &wk, &ap[kk]);
	wk.r = z__1.r, wk.i = z__1.i;
	z_div(&z__1, &wkm, &ap[kk]);
	wkm.r = z__1.r, wkm.i = z__1.i;
	kp1 = k + 1;
	kj = kk + k;
	if (kp1 > *n) {
	    goto L100;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    d_cnjg(&z__4, &ap[kj]);
	    z__3.r = wkm.r * z__4.r - wkm.i * z__4.i, z__3.i = wkm.r * z__4.i 
		    + wkm.i * z__4.r;
	    z__2.r = z__[i__3].r + z__3.r, z__2.i = z__[i__3].i + z__3.i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    sm += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(
		    d__2));
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &ap[kj]);
	    z__2.r = wk.r * z__3.r - wk.i * z__3.i, z__2.i = wk.r * z__3.i + 
		    wk.i * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = j;
	    s += (d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&z__[j]), 
		    abs(d__2));
	    kj += j;
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	z__1.r = wkm.r - wk.r, z__1.i = wkm.i - wk.i;
	t.r = z__1.r, t.i = z__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	kj = kk + k;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &ap[kj]);
	    z__2.r = t.r * z__3.r - t.i * z__3.i, z__2.i = t.r * z__3.i + t.i 
		    * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    kj += j;
/* L80: */
	}
L90:
L100:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L110: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*        SOLVE R*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = kk;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= ap[i__3].r) {
	    goto L120;
	}
	i__2 = kk;
	i__3 = k;
	s = ap[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		z__[k]), abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
L120:
	i__2 = k;
	z_div(&z__1, &z__[k], &ap[kk]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	kk -= k;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
/* L130: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*        SOLVE CTRANS(R)*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	i__3 = k;
	i__4 = k - 1;
	cdotc_(&z__2, &i__4, &ap[kk + 1], &c__1, &z__[1], &c__1);
	z__1.r = z__[i__3].r - z__2.r, z__1.i = z__[i__3].i - z__2.i;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	kk += k;
	i__2 = k;
	i__3 = kk;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= ap[i__3].r) {
	    goto L140;
	}
	i__2 = kk;
	i__3 = k;
	s = ap[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		z__[k]), abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	i__2 = k;
	z_div(&z__1, &z__[k], &ap[kk]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
/* L150: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        SOLVE R*Z = V */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = kk;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= ap[i__3].r) {
	    goto L160;
	}
	i__2 = kk;
	i__3 = k;
	s = ap[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		z__[k]), abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	i__2 = k;
	z_div(&z__1, &z__[k], &ap[kk]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	kk -= k;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
/* L170: */
    }
/*        MAKE ZNORM = 1.0 */
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
L180:
    return 0;
} /* cppco_ */

/* Subroutine */ int cppfa_(doublecomplex *ap, integer *n, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal s;
    static doublecomplex t;
    static integer jj, kj, kk, jm1;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CPPFA FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX */
/*     STORED IN PACKED FORM. */

/*     CPPFA IS USUALLY CALLED BY CPPCO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */
/*     (TIME FOR CPPCO) = (1 + 18/N)*(TIME FOR CPPFA) . */

/*     ON ENTRY */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                THE PACKED FORM OF A HERMITIAN MATRIX  A .  THE */
/*                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY */
/*                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 . */
/*                SEE COMMENTS BELOW FOR DETAILS. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        AP      AN UPPER TRIANGULAR MATRIX  R , STORED IN PACKED */
/*                FORM, SO THAT  A = CTRANS(R)*R . */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT */
/*                     POSITIVE DEFINITE. */


/*     PACKED STORAGE */

/*          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER */
/*          TRIANGLE OF A HERMITIAN MATRIX. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K) = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CDOTC */
/*     FORTRAN DIMAG,DCMPLX,CONJG,DBLE,SQRT */

/*     INTERNAL VARIABLES */

/*     BEGIN BLOCK WITH ...EXITS TO 40 */


    /* Parameter adjustments */
    --ap;

    /* Function Body */
    jj = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	kj = jj;
	kk = 0;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    ++kj;
	    i__3 = kj;
	    i__4 = k - 1;
	    cdotc_(&z__2, &i__4, &ap[kk + 1], &c__1, &ap[jj + 1], &c__1);
	    z__1.r = ap[i__3].r - z__2.r, z__1.i = ap[i__3].i - z__2.i;
	    t.r = z__1.r, t.i = z__1.i;
	    kk += k;
	    z_div(&z__1, &t, &ap[kk]);
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = kj;
	    ap[i__3].r = t.r, ap[i__3].i = t.i;
	    d_cnjg(&z__2, &t);
	    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i 
		    * z__2.r;
	    s += z__1.r;
/* L10: */
	}
L20:
	jj += j;
	i__2 = jj;
	s = ap[i__2].r - s;
/*     ......EXIT */
	if (s <= 0. || d_imag(&ap[jj]) != 0.) {
	    goto L40;
	}
	i__2 = jj;
	d__1 = sqrt(s);
	z__1.r = d__1, z__1.i = 0.;
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* cppfa_ */

/* Subroutine */ int cppsl_(doublecomplex *ap, integer *n, doublecomplex *b)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex t;
    static integer kb, kk;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CPPSL SOLVES THE COMPLEX HERMITIAN POSITIVE DEFINITE SYSTEM */
/*     A * X = B */
/*     USING THE FACTORS COMPUTED BY CPPCO OR CPPFA. */

/*     ON ENTRY */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                THE OUTPUT FROM CPPCO OR CPPFA. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        B       COMPLEX(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES */
/*        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE */
/*        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED */
/*        CORRECTLY AND  INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL CPPCO(AP,N,RCOND,Z,INFO) */
/*           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CPPSL(AP,N,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CDOTC */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    --b;
    --ap;

    /* Function Body */
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	cdotc_(&z__1, &i__2, &ap[kk + 1], &c__1, &b[1], &c__1);
	t.r = z__1.r, t.i = z__1.i;
	kk += k;
	i__2 = k;
	i__3 = k;
	z__2.r = b[i__3].r - t.r, z__2.i = b[i__3].i - t.i;
	z_div(&z__1, &z__2, &ap[kk]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L10: */
    }
    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	z_div(&z__1, &b[k], &ap[kk]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	kk -= k;
	i__2 = k;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &ap[kk + 1], &c__1, &b[1], &c__1);
/* L20: */
    }
    return 0;
} /* cppsl_ */

/* Subroutine */ int cppdi_(doublecomplex *ap, integer *n, doublereal *det, 
	integer *job)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static doublecomplex t;
    static integer j1, k1, ii, jj, kj, kk, jm1, kp1;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CPPDI COMPUTES THE DETERMINANT AND INVERSE */
/*     OF A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX */
/*     USING THE FACTORS COMPUTED BY CPPCO OR CPPFA . */

/*     ON ENTRY */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                THE OUTPUT FROM CPPCO OR CPPFA. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        JOB     INTEGER */
/*                = 11   BOTH DETERMINANT AND INVERSE. */
/*                = 01   INVERSE ONLY. */
/*                = 10   DETERMINANT ONLY. */

/*     ON RETURN */

/*        AP      THE UPPER TRIANGULAR HALF OF THE INVERSE . */
/*                THE STRICT LOWER TRIANGLE IS UNALTERED. */

/*        DET     REAL(2) */
/*                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED. */
/*                OTHERWISE NOT REFERENCED. */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. DET(1) .LT. 10.0 */
/*                OR  DET(1) .EQ. 0.0 . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL AND THE INVERSE IS REQUESTED. */
/*        IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED CORRECTLY */
/*        AND IF CPOCO OR CPOFA HAS SET INFO .EQ. 0 . */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CSCAL */
/*     FORTRAN CONJG,MOD,DBLE */

/*     INTERNAL VARIABLES */


/*     COMPUTE DETERMINANT */

    /* Parameter adjustments */
    --det;
    --ap;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    s = 10.;
    ii = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii += i__;
	i__2 = ii;
/* Computing 2nd power */
	d__1 = ap[i__2].r;
	det[1] = d__1 * d__1 * det[1];
/*        ...EXIT */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (det[1] >= 1.) {
	    goto L20;
	}
	det[1] = s * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (det[1] < s) {
	    goto L40;
	}
	det[1] /= s;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*     COMPUTE INVERSE(R) */

    if (*job % 10 == 0) {
	goto L140;
    }
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	k1 = kk + 1;
	kk += k;
	i__2 = kk;
	z_div(&z__1, &c_b1092, &ap[kk]);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = kk;
	z__1.r = -ap[i__2].r, z__1.i = -ap[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	cscal_(&i__2, &t, &ap[k1], &c__1);
	kp1 = k + 1;
	j1 = kk + 1;
	kj = kk + k;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = kj;
	    t.r = ap[i__3].r, t.i = ap[i__3].i;
	    i__3 = kj;
	    ap[i__3].r = 0., ap[i__3].i = 0.;
	    caxpy_(&k, &t, &ap[k1], &c__1, &ap[j1], &c__1);
	    j1 += j;
	    kj += j;
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        FORM  INVERSE(R) * CTRANS(INVERSE(R)) */

    jj = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	j1 = jj + 1;
	jj += j;
	jm1 = j - 1;
	k1 = 1;
	kj = j1;
	if (jm1 < 1) {
	    goto L120;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    d_cnjg(&z__1, &ap[kj]);
	    t.r = z__1.r, t.i = z__1.i;
	    caxpy_(&k, &t, &ap[j1], &c__1, &ap[k1], &c__1);
	    k1 += k;
	    ++kj;
/* L110: */
	}
L120:
	d_cnjg(&z__1, &ap[jj]);
	t.r = z__1.r, t.i = z__1.i;
	cscal_(&j, &t, &ap[j1], &c__1);
/* L130: */
    }
L140:
    return 0;
} /* cppdi_ */

/* Subroutine */ int cpbco_(doublecomplex *abd, integer *lda, integer *n, 
	integer *m, doublereal *rcond, doublecomplex *z__, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal s;
    static doublecomplex t;
    static integer j2, kb, la, lb;
    static doublecomplex ek;
    static integer lm;
    static doublereal sm;
    static doublecomplex wk;
    static integer mu, kp1;
    static doublecomplex wkm;
    extern /* Subroutine */ int cpbfa_(doublecomplex *, integer *, integer *, 
	    integer *, integer *);
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal anorm;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);


/*     CPBCO FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX */
/*     STORED IN BAND FORM AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, CPBFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW CPBCO BY CPBSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CPBCO BY CPBSL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CPBCO BY CPBDI. */

/*     ON ENTRY */

/*        ABD     COMPLEX(LDA, N) */
/*                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER */
/*                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE */
/*                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE */
/*                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */
/*                LDA MUST BE .GE. M + 1 . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        M       INTEGER */
/*                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */
/*                0 .LE. M .LT. N . */

/*     ON RETURN */

/*        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND */
/*                FORM, SO THAT  A = CTRANS(R)*R . */
/*                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE. */

/*        RCOND   REAL */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS.  IF INFO .NE. 0 , RCOND IS UNCHANGED. */

/*        Z       COMPLEX(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS SINGULAR TO WORKING PRECISION, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */
/*                IF  INFO .NE. 0 , Z  IS UNCHANGED. */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR */
/*                     OF ORDER  K  IS NOT POSITIVE DEFINITE. */

/*     BAND STORAGE */

/*           IF  A  IS A HERMITIAN POSITIVE DEFINITE BAND MATRIX, */
/*           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT. */

/*                   M = (BAND WIDTH ABOVE DIAGONAL) */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX0(1, J-M) */
/*                      DO 10 I = I1, J */
/*                         K = I-J+M+1 */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*           THIS USES  M + 1  ROWS OF  A , EXCEPT FOR THE  M BY M */
/*           UPPER LEFT TRIANGLE, WHICH IS IGNORED. */

/*     EXAMPLE..  IF THE ORIGINAL MATRIX IS */

/*           11 12 13  0  0  0 */
/*           12 22 23 24  0  0 */
/*           13 23 33 34 35  0 */
/*            0 24 34 44 45 46 */
/*            0  0 35 45 55 56 */
/*            0  0  0 46 56 66 */

/*     THEN  N = 6 , M = 2  AND  ABD  SHOULD CONTAIN */

/*            *  * 13 24 35 46 */
/*            * 12 23 34 45 56 */
/*           11 22 33 44 55 66 */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK CPBFA */
/*     BLAS CAXPY,CDOTC,CSSCAL,SCASUM */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,CONJG,MAX0,MIN0,DBLE */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	i__2 = j, i__3 = *m + 1;
	l = min(i__2,i__3);
/* Computing MAX */
	i__2 = *m + 2 - j;
	mu = max(i__2,1);
	i__2 = j;
	d__1 = scasum_(&l, &abd[mu + j * abd_dim1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	k = j - l;
	if (*m < mu) {
	    goto L20;
	}
	i__2 = *m;
	for (i__ = mu; i__ <= i__2; ++i__) {
	    ++k;
	    i__3 = k;
	    i__4 = k;
	    i__5 = i__ + j * abd_dim1;
	    d__3 = z__[i__4].r + ((d__1 = abd[i__5].r, abs(d__1)) + (d__2 = 
		    d_imag(&abd[i__ + j * abd_dim1]), abs(d__2)));
	    z__1.r = d__3, z__1.i = 0.;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = j;
	d__1 = anorm, d__2 = z__[i__2].r;
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    cpbfa_(&abd[abd_offset], lda, n, m, info);
    if (*info != 0) {
	goto L180;
    }

/*        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*        GROWTH IN THE ELEMENTS OF W  WHERE  CTRANS(R)*W = E . */
/*        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*        SOLVE CTRANS(R)*W = E */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L50: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) != 0.) {
	    i__3 = k;
	    z__2.r = -z__[i__3].r, z__2.i = -z__[i__3].i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	    d__8 = (d__5 = z__1.r, abs(d__5)) + (d__6 = d_imag(&z__1), abs(
		    d__6));
	    z__4.r = z__1.r / d__8, z__4.i = z__1.i / d__8;
	    z__3.r = d__7 * z__4.r, z__3.i = d__7 * z__4.i;
	    ek.r = z__3.r, ek.i = z__3.i;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = *m + 1 + k * abd_dim1;
	if ((d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(d__2)) <= 
		abd[i__3].r) {
	    goto L60;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = *m + 1 + k * abd_dim1;
	s = abd[i__3].r / ((d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1),
		 abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
	z__2.r = s, z__2.i = 0.;
	z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + 
		z__2.i * ek.r;
	ek.r = z__1.r, ek.i = z__1.i;
L60:
	i__2 = k;
	z__1.r = ek.r - z__[i__2].r, z__1.i = ek.i - z__[i__2].i;
	wk.r = z__1.r, wk.i = z__1.i;
	z__2.r = -ek.r, z__2.i = -ek.i;
	i__2 = k;
	z__1.r = z__2.r - z__[i__2].r, z__1.i = z__2.i - z__[i__2].i;
	wkm.r = z__1.r, wkm.i = z__1.i;
	s = (d__1 = wk.r, abs(d__1)) + (d__2 = d_imag(&wk), abs(d__2));
	sm = (d__1 = wkm.r, abs(d__1)) + (d__2 = d_imag(&wkm), abs(d__2));
	z_div(&z__1, &wk, &abd[*m + 1 + k * abd_dim1]);
	wk.r = z__1.r, wk.i = z__1.i;
	z_div(&z__1, &wkm, &abd[*m + 1 + k * abd_dim1]);
	wkm.r = z__1.r, wkm.i = z__1.i;
	kp1 = k + 1;
/* Computing MIN */
	i__2 = k + *m;
	j2 = min(i__2,*n);
	i__ = *m + 1;
	if (kp1 > j2) {
	    goto L100;
	}
	i__2 = j2;
	for (j = kp1; j <= i__2; ++j) {
	    --i__;
	    i__3 = j;
	    d_cnjg(&z__4, &abd[i__ + j * abd_dim1]);
	    z__3.r = wkm.r * z__4.r - wkm.i * z__4.i, z__3.i = wkm.r * z__4.i 
		    + wkm.i * z__4.r;
	    z__2.r = z__[i__3].r + z__3.r, z__2.i = z__[i__3].i + z__3.i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    sm += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(
		    d__2));
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &abd[i__ + j * abd_dim1]);
	    z__2.r = wk.r * z__3.r - wk.i * z__3.i, z__2.i = wk.r * z__3.i + 
		    wk.i * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = j;
	    s += (d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&z__[j]), 
		    abs(d__2));
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	z__1.r = wkm.r - wk.r, z__1.i = wkm.i - wk.i;
	t.r = z__1.r, t.i = z__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	i__ = *m + 1;
	i__2 = j2;
	for (j = kp1; j <= i__2; ++j) {
	    --i__;
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &abd[i__ + j * abd_dim1]);
	    z__2.r = t.r * z__3.r - t.i * z__3.i, z__2.i = t.r * z__3.i + t.i 
		    * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
/* L80: */
	}
L90:
L100:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L110: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*        SOLVE  R*Y = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = *m + 1 + k * abd_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= abd[i__3].r) {
	    goto L120;
	}
	i__2 = *m + 1 + k * abd_dim1;
	i__3 = k;
	s = abd[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		z__[k]), abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
L120:
	i__2 = k;
	z_div(&z__1, &z__[k], &abd[*m + 1 + k * abd_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	caxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L130: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*        SOLVE CTRANS(R)*V = Y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	i__2 = k;
	i__3 = k;
	cdotc_(&z__2, &lm, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
	z__1.r = z__[i__3].r - z__2.r, z__1.i = z__[i__3].i - z__2.i;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	i__2 = k;
	i__3 = *m + 1 + k * abd_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= abd[i__3].r) {
	    goto L140;
	}
	i__2 = *m + 1 + k * abd_dim1;
	i__3 = k;
	s = abd[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		z__[k]), abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	i__2 = k;
	z_div(&z__1, &z__[k], &abd[*m + 1 + k * abd_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
/* L150: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        SOLVE  R*Z = W */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = *m + 1 + k * abd_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= abd[i__3].r) {
	    goto L160;
	}
	i__2 = *m + 1 + k * abd_dim1;
	i__3 = k;
	s = abd[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&
		z__[k]), abs(d__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	i__2 = k;
	z_div(&z__1, &z__[k], &abd[*m + 1 + k * abd_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	caxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L170: */
    }
/*        MAKE ZNORM = 1.0 */
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
L180:
    return 0;
} /* cpbco_ */

/* Subroutine */ int cpbfa_(doublecomplex *abd, integer *lda, integer *n, 
	integer *m, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);
    double d_imag(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal s;
    static doublecomplex t;
    static integer ik, jk, mu;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CPBFA FACTORS A COMPLEX HERMITIAN POSITIVE DEFINITE MATRIX */
/*     STORED IN BAND FORM. */

/*     CPBFA IS USUALLY CALLED BY CPBCO, BUT IT CAN BE CALLED */
/*     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED. */

/*     ON ENTRY */

/*        ABD     COMPLEX(LDA, N) */
/*                THE MATRIX TO BE FACTORED.  THE COLUMNS OF THE UPPER */
/*                TRIANGLE ARE STORED IN THE COLUMNS OF ABD AND THE */
/*                DIAGONALS OF THE UPPER TRIANGLE ARE STORED IN THE */
/*                ROWS OF ABD .  SEE THE COMMENTS BELOW FOR DETAILS. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */
/*                LDA MUST BE .GE. M + 1 . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        M       INTEGER */
/*                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */
/*                0 .LE. M .LT. N . */

/*     ON RETURN */

/*        ABD     AN UPPER TRIANGULAR MATRIX  R , STORED IN BAND */
/*                FORM, SO THAT  A = CTRANS(R)*R . */

/*        INFO    INTEGER */
/*                = 0  FOR NORMAL RETURN. */
/*                = K  IF THE LEADING MINOR OF ORDER  K  IS NOT */
/*                     POSITIVE DEFINITE. */

/*     BAND STORAGE */

/*           IF  A  IS A HERMITIAN POSITIVE DEFINITE BAND MATRIX, */
/*           THE FOLLOWING PROGRAM SEGMENT WILL SET UP THE INPUT. */

/*                   M = (BAND WIDTH ABOVE DIAGONAL) */
/*                   DO 20 J = 1, N */
/*                      I1 = MAX0(1, J-M) */
/*                      DO 10 I = I1, J */
/*                         K = I-J+M+1 */
/*                         ABD(K,J) = A(I,J) */
/*                10    CONTINUE */
/*                20 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CDOTC */
/*     FORTRAN DIMAG,DCMPLX,CONJG,MAX0,DBLE,SQRT */

/*     INTERNAL VARIABLES */

/*     BEGIN BLOCK WITH ...EXITS TO 40 */


    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	ik = *m + 1;
/* Computing MAX */
	i__2 = j - *m;
	jk = max(i__2,1);
/* Computing MAX */
	i__2 = *m + 2 - j;
	mu = max(i__2,1);
	if (*m < mu) {
	    goto L20;
	}
	i__2 = *m;
	for (k = mu; k <= i__2; ++k) {
	    i__3 = k + j * abd_dim1;
	    i__4 = k - mu;
	    cdotc_(&z__2, &i__4, &abd[ik + jk * abd_dim1], &c__1, &abd[mu + j 
		    * abd_dim1], &c__1);
	    z__1.r = abd[i__3].r - z__2.r, z__1.i = abd[i__3].i - z__2.i;
	    t.r = z__1.r, t.i = z__1.i;
	    z_div(&z__1, &t, &abd[*m + 1 + jk * abd_dim1]);
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = k + j * abd_dim1;
	    abd[i__3].r = t.r, abd[i__3].i = t.i;
	    d_cnjg(&z__2, &t);
	    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i 
		    * z__2.r;
	    s += z__1.r;
	    --ik;
	    ++jk;
/* L10: */
	}
L20:
	i__2 = *m + 1 + j * abd_dim1;
	s = abd[i__2].r - s;
/*     ......EXIT */
	if (s <= 0. || d_imag(&abd[*m + 1 + j * abd_dim1]) != 0.) {
	    goto L40;
	}
	i__2 = *m + 1 + j * abd_dim1;
	d__1 = sqrt(s);
	z__1.r = d__1, z__1.i = 0.;
	abd[i__2].r = z__1.r, abd[i__2].i = z__1.i;
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* cpbfa_ */

/* Subroutine */ int cpbsl_(doublecomplex *abd, integer *lda, integer *n, 
	integer *m, doublecomplex *b)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex t;
    static integer kb, la, lb, lm;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CPBSL SOLVES THE COMPLEX HERMITIAN POSITIVE DEFINITE BAND */
/*     SYSTEM  A*X = B */
/*     USING THE FACTORS COMPUTED BY CPBCO OR CPBFA. */

/*     ON ENTRY */

/*        ABD     COMPLEX(LDA, N) */
/*                THE OUTPUT FROM CPBCO OR CPBFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        M       INTEGER */
/*                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */

/*        B       COMPLEX(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS */
/*        A ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES */
/*        SINGULARITY BUT IT IS USUALLY CAUSED BY IMPROPER SUBROUTINE */
/*        ARGUMENTS.  IT WILL NOT OCCUR IF THE SUBROUTINES ARE CALLED */
/*        CORRECTLY AND  INFO .EQ. 0 . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL CPBCO(ABD,LDA,N,RCOND,Z,INFO) */
/*           IF (RCOND IS TOO SMALL .OR. INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CPBSL(ABD,LDA,N,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CDOTC */
/*     FORTRAN MIN0 */

/*     INTERNAL VARIABLES */


/*     SOLVE CTRANS(R)*Y = B */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	cdotc_(&z__1, &lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k;
	i__3 = k;
	z__2.r = b[i__3].r - t.r, z__2.i = b[i__3].i - t.i;
	z_div(&z__1, &z__2, &abd[*m + 1 + k * abd_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L10: */
    }

/*     SOLVE R*X = Y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	i__2 = k;
	z_div(&z__1, &b[k], &abd[*m + 1 + k * abd_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	i__2 = k;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	caxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/* L20: */
    }
    return 0;
} /* cpbsl_ */

/* Subroutine */ int cpbdi_(doublecomplex *abd, integer *lda, integer *n, 
	integer *m, doublereal *det)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal s;


/*     CPBDI COMPUTES THE DETERMINANT */
/*     OF A COMPLEX HERMITIAN POSITIVE DEFINITE BAND MATRIX */
/*     USING THE FACTORS COMPUTED BY CPBCO OR CPBFA. */
/*     IF THE INVERSE IS NEEDED, USE CPBSL  N  TIMES. */

/*     ON ENTRY */

/*        ABD     COMPLEX(LDA, N) */
/*                THE OUTPUT FROM CPBCO OR CPBFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  ABD . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        M       INTEGER */
/*                THE NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL. */

/*     ON RETURN */

/*        DET     REAL(2) */
/*                DETERMINANT OF ORIGINAL MATRIX IN THE FORM */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. DET(1) .LT. 10.0 */
/*                OR  DET(1) .EQ. 0.0 . */

/*     LINPACK.  THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     FORTRAN DBLE */

/*     INTERNAL VARIABLES */


/*     COMPUTE DETERMINANT */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --det;

    /* Function Body */
    det[1] = 1.;
    det[2] = 0.;
    s = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m + 1 + i__ * abd_dim1;
/* Computing 2nd power */
	d__1 = abd[i__2].r;
	det[1] = d__1 * d__1 * det[1];
/*     ...EXIT */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (det[1] >= 1.) {
	    goto L20;
	}
	det[1] = s * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (det[1] < s) {
	    goto L40;
	}
	det[1] /= s;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
    return 0;
} /* cpbdi_ */

/* Subroutine */ int csico_(doublecomplex *a, integer *lda, integer *n, 
	integer *kpvt, doublereal *rcond, doublecomplex *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static doublecomplex t, ak, bk, ek;
    static integer kp, ks, jm1, kps;
    static doublecomplex akm1, bkm1;
    static integer info;
    extern /* Subroutine */ int csifa_(doublecomplex *, integer *, integer *, 
	    integer *, integer *);
    static doublecomplex denom;
    static doublereal anorm;
    extern /* Double Complex */ VOID cdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);


/*     CSICO FACTORS A COMPLEX SYMMETRIC MATRIX BY ELIMINATION WITH */
/*     SYMMETRIC PIVOTING AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, CSIFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW CSICO BY CSISL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CSICO BY CSISL. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW CSICO BY CSIDI. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CSICO BY CSIDI. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA, N) */
/*                THE SYMMETRIC MATRIX TO BE FACTORED. */
/*                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE */
/*                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        RCOND   REAL */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       COMPLEX(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK CSIFA */
/*     BLAS CAXPY,CDOTU,CSSCAL,SCASUM */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,IABS,DBLE */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A USING ONLY UPPER HALF */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	d__1 = scasum_(&j, &a[j * a_dim1 + 1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__ + j * a_dim1;
	    d__3 = z__[i__4].r + ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[i__ + j * a_dim1]), abs(d__2)));
	    z__1.r = d__3, z__1.i = 0.;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = j;
	d__1 = anorm, d__2 = z__[i__2].r;
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    csifa_(&a[a_offset], lda, n, &kpvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE U*D*W = E */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L50: */
    }
    k = *n;
L60:
    if (k == 0) {
	goto L120;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L70;
    }
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L70:
    i__1 = k;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    != 0.) {
	d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	i__2 = k;
	i__3 = k;
	d__8 = (d__5 = z__[i__3].r, abs(d__5)) + (d__6 = d_imag(&z__[k]), abs(
		d__6));
	z__2.r = z__[i__2].r / d__8, z__2.i = z__[i__2].i / d__8;
	z__1.r = d__7 * z__2.r, z__1.i = d__7 * z__2.i;
	ek.r = z__1.r, ek.i = z__1.i;
    }
    i__1 = k;
    i__2 = k;
    z__1.r = z__[i__2].r + ek.r, z__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 1) {
	goto L80;
    }
    i__1 = k - 1;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k - 1]), abs(
	    d__2)) != 0.) {
	d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	i__2 = k - 1;
	i__3 = k - 1;
	d__8 = (d__5 = z__[i__3].r, abs(d__5)) + (d__6 = d_imag(&z__[k - 1]), 
		abs(d__6));
	z__2.r = z__[i__2].r / d__8, z__2.i = z__[i__2].i / d__8;
	z__1.r = d__7 * z__2.r, z__1.i = d__7 * z__2.i;
	ek.r = z__1.r, ek.i = z__1.i;
    }
    i__1 = k - 1;
    i__2 = k - 1;
    z__1.r = z__[i__2].r + ek.r, z__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
	    c__1);
L80:
    if (ks == 2) {
	goto L100;
    }
    i__1 = k;
    i__2 = k + k * a_dim1;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    <= (d__3 = a[i__2].r, abs(d__3)) + (d__4 = d_imag(&a[k + k * 
	    a_dim1]), abs(d__4))) {
	goto L90;
    }
    i__1 = k + k * a_dim1;
    i__2 = k;
    s = ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2))) / ((d__3 = z__[i__2].r, abs(d__3)) + (d__4 = d_imag(&
	    z__[k]), abs(d__4)));
    csscal_(n, &s, &z__[1], &c__1);
    z__2.r = s, z__2.i = 0.;
    z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + z__2.i * 
	    ek.r;
    ek.r = z__1.r, ek.i = z__1.i;
L90:
    i__1 = k + k * a_dim1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2)) != 0.) {
	i__2 = k;
	z_div(&z__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
    }
    i__1 = k + k * a_dim1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2)) == 0.) {
	i__2 = k;
	z__[i__2].r = 1., z__[i__2].i = 0.;
    }
    goto L110;
L100:
    z_div(&z__1, &a[k + k * a_dim1], &a[k - 1 + k * a_dim1]);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &a[k - 1 + k * a_dim1]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    z_div(&z__1, &z__[k], &a[k - 1 + k * a_dim1]);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &z__[k - 1], &a[k - 1 + k * a_dim1]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
L110:
    k -= ks;
    goto L60;
L120:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*     SOLVE TRANS(U)*Y = W */

    k = 1;
L130:
    if (k > *n) {
	goto L160;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L150;
    }
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotu_(&z__2, &i__3, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &c__1);
	z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
	z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L140;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L140:
L150:
    k += ks;
    goto L130;
L160:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE U*D*V = Y */

    k = *n;
L170:
    if (k == 0) {
	goto L230;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == ks) {
	goto L190;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L180;
    }
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L180:
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - ks;
	caxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
L190:
    if (ks == 2) {
	goto L210;
    }
    i__1 = k;
    i__2 = k + k * a_dim1;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    <= (d__3 = a[i__2].r, abs(d__3)) + (d__4 = d_imag(&a[k + k * 
	    a_dim1]), abs(d__4))) {
	goto L200;
    }
    i__1 = k + k * a_dim1;
    i__2 = k;
    s = ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2))) / ((d__3 = z__[i__2].r, abs(d__3)) + (d__4 = d_imag(&
	    z__[k]), abs(d__4)));
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;
L200:
    i__1 = k + k * a_dim1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2)) != 0.) {
	i__2 = k;
	z_div(&z__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
    }
    i__1 = k + k * a_dim1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2)) == 0.) {
	i__2 = k;
	z__[i__2].r = 1., z__[i__2].i = 0.;
    }
    goto L220;
L210:
    z_div(&z__1, &a[k + k * a_dim1], &a[k - 1 + k * a_dim1]);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &a[k - 1 + k * a_dim1]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    z_div(&z__1, &z__[k], &a[k - 1 + k * a_dim1]);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &z__[k - 1], &a[k - 1 + k * a_dim1]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
L220:
    k -= ks;
    goto L170;
L230:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE TRANS(U)*Z = V */

    k = 1;
L240:
    if (k > *n) {
	goto L270;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L260;
    }
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotu_(&z__2, &i__3, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &c__1);
	z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
	z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L250;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L250:
L260:
    k += ks;
    goto L240;
L270:
/*     MAKE ZNORM = 1.0 */
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* csico_ */

/* Subroutine */ int csifa_(doublecomplex *a, integer *lda, integer *n, 
	integer *kpvt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublecomplex t, ak, bk;
    static integer jj, km1, km2;
    static doublecomplex akm1, bkm1;
    static integer imax, jmax;
    static doublecomplex mulk;
    static logical swap;
    static doublereal alpha;
    static doublecomplex denom;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer kstep, imaxp1;
    static doublecomplex mulkm1;
    static doublereal absakk;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    static doublereal colmax, rowmax;


/*     CSIFA FACTORS A COMPLEX SYMMETRIC MATRIX BY ELIMINATION */
/*     WITH SYMMETRIC PIVOTING. */

/*     TO SOLVE  A*X = B , FOLLOW CSIFA BY CSISL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CSIFA BY CSISL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CSIFA BY CSIDI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW CSIFA BY CSIDI. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA,N) */
/*                THE SYMMETRIC MATRIX TO BE FACTORED. */
/*                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE */
/*                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS */
/*                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE, */
/*                     BUT IT DOES INDICATE THAT CSISL OR CSIDI MAY */
/*                     DIVIDE BY ZERO IF CALLED. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CSWAP,ICAMAX */
/*     FORTRAN ABS,DIMAG,MAX,DBLE,SQRT */

/*     INTERNAL VARIABLES */



/*     INITIALIZE */

/*     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;

    /* Function Body */
    alpha = (sqrt(17.) + 1.) / 8.;

    *info = 0;

/*     MAIN LOOP ON K, WHICH GOES FROM N TO 1. */

    k = *n;
L10:

/*        LEAVE THE LOOP IF K=0 OR K=1. */

/*     ...EXIT */
    if (k == 0) {
	goto L200;
    }
    if (k > 1) {
	goto L20;
    }
    kpvt[1] = 1;
    i__1 = a_dim1 + 1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[a_dim1 + 1]), abs(
	    d__2)) == 0.) {
	*info = 1;
    }
/*     ......EXIT */
    goto L200;
L20:

/*        THIS SECTION OF CODE DETERMINES THE KIND OF */
/*        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED, */
/*        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND */
/*        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS */
/*        REQUIRED. */

    km1 = k - 1;
    i__1 = k + k * a_dim1;
    absakk = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]
	    ), abs(d__2));

/*        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*        COLUMN K. */

    i__1 = k - 1;
    imax = icamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = imax + k * a_dim1;
    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + k * 
	    a_dim1]), abs(d__2));
    if (absakk < alpha * colmax) {
	goto L30;
    }
    kstep = 1;
    swap = FALSE_;
    goto L90;
L30:

/*           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*           ROW IMAX. */

    rowmax = 0.;
    imaxp1 = imax + 1;
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = imax + j * a_dim1;
	d__3 = rowmax, d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&
		a[imax + j * a_dim1]), abs(d__2));
	rowmax = max(d__3,d__4);
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = icamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
    i__1 = jmax + imax * a_dim1;
    d__3 = rowmax, d__4 = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
	    jmax + imax * a_dim1]), abs(d__2));
    rowmax = max(d__3,d__4);
L50:
    i__1 = imax + imax * a_dim1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + imax * 
	    a_dim1]), abs(d__2)) < alpha * rowmax) {
	goto L60;
    }
    kstep = 1;
    swap = TRUE_;
    goto L80;
L60:
    if (absakk < alpha * colmax * (colmax / rowmax)) {
	goto L70;
    }
    kstep = 1;
    swap = FALSE_;
    goto L80;
L70:
    kstep = 2;
    swap = imax != km1;
L80:
L90:
    if (max(absakk,colmax) != 0.) {
	goto L100;
    }

/*           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP. */

    kpvt[k] = k;
    *info = k;
    goto L190;
L100:
    if (kstep == 2) {
	goto L140;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (! swap) {
	goto L120;
    }

/*              PERFORM AN INTERCHANGE. */

    cswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	i__2 = j + k * a_dim1;
	t.r = a[i__2].r, t.i = a[i__2].i;
	i__2 = j + k * a_dim1;
	i__3 = imax + j * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = imax + j * a_dim1;
	a[i__2].r = t.r, a[i__2].i = t.i;
/* L110: */
    }
L120:

/*           PERFORM THE ELIMINATION. */

    i__1 = km1;
    for (jj = 1; jj <= i__1; ++jj) {
	j = k - jj;
	i__2 = j + k * a_dim1;
	z__2.r = -a[i__2].r, z__2.i = -a[i__2].i;
	z_div(&z__1, &z__2, &a[k + k * a_dim1]);
	mulk.r = z__1.r, mulk.i = z__1.i;
	t.r = mulk.r, t.i = mulk.i;
	caxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	i__2 = j + k * a_dim1;
	a[i__2].r = mulk.r, a[i__2].i = mulk.i;
/* L130: */
    }

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = k;
    if (swap) {
	kpvt[k] = imax;
    }
    goto L190;
L140:

/*           2 X 2 PIVOT BLOCK. */

    if (! swap) {
	goto L160;
    }

/*              PERFORM AN INTERCHANGE. */

    cswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[(k - 1) * a_dim1 + 1], &
	    c__1);
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	i__2 = j + (k - 1) * a_dim1;
	t.r = a[i__2].r, t.i = a[i__2].i;
	i__2 = j + (k - 1) * a_dim1;
	i__3 = imax + j * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = imax + j * a_dim1;
	a[i__2].r = t.r, a[i__2].i = t.i;
/* L150: */
    }
    i__1 = k - 1 + k * a_dim1;
    t.r = a[i__1].r, t.i = a[i__1].i;
    i__1 = k - 1 + k * a_dim1;
    i__2 = imax + k * a_dim1;
    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
    i__1 = imax + k * a_dim1;
    a[i__1].r = t.r, a[i__1].i = t.i;
L160:

/*           PERFORM THE ELIMINATION. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    z_div(&z__1, &a[k + k * a_dim1], &a[k - 1 + k * a_dim1]);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &a[k - 1 + k * a_dim1]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = 1. - z__2.r, z__1.i = -z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	z_div(&z__1, &a[j + k * a_dim1], &a[k - 1 + k * a_dim1]);
	bk.r = z__1.r, bk.i = z__1.i;
	z_div(&z__1, &a[j + (k - 1) * a_dim1], &a[k - 1 + k * a_dim1]);
	bkm1.r = z__1.r, bkm1.i = z__1.i;
	z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + 
		akm1.i * bk.r;
	z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
	z_div(&z__1, &z__2, &denom);
	mulk.r = z__1.r, mulk.i = z__1.i;
	z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i 
		* bkm1.r;
	z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
	z_div(&z__1, &z__2, &denom);
	mulkm1.r = z__1.r, mulkm1.i = z__1.i;
	t.r = mulk.r, t.i = mulk.i;
	caxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	t.r = mulkm1.r, t.i = mulkm1.i;
	caxpy_(&j, &t, &a[(k - 1) * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		c__1);
	i__2 = j + k * a_dim1;
	a[i__2].r = mulk.r, a[i__2].i = mulk.i;
	i__2 = j + (k - 1) * a_dim1;
	a[i__2].r = mulkm1.r, a[i__2].i = mulkm1.i;
/* L170: */
    }
L180:

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = 1 - k;
    if (swap) {
	kpvt[k] = -imax;
    }
    kpvt[k - 1] = kpvt[k];
L190:
    k -= kstep;
    goto L10;
L200:
    return 0;
} /* csifa_ */

/* Subroutine */ int csisl_(doublecomplex *a, integer *lda, integer *n, 
	integer *kpvt, doublecomplex *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex ak, bk;
    static integer kp;
    static doublecomplex akm1, bkm1, temp, denom;
    extern /* Double Complex */ VOID cdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CSISL SOLVES THE COMPLEX SYMMETRIC SYSTEM */
/*     A * X = B */
/*     USING THE FACTORS COMPUTED BY CSIFA. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA,N) */
/*                THE OUTPUT FROM CSIFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CSIFA. */

/*        B       COMPLEX(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO MAY OCCUR IF  CSICO  HAS SET RCOND .EQ. 0.0 */
/*        OR  CSIFA  HAS SET INFO .NE. 0  . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL CSIFA(A,LDA,N,KPVT,INFO) */
/*           IF (INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CSISL(A,LDA,N,KPVT,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CDOTU */
/*     FORTRAN IABS */

/*     INTERNAL VARIABLES. */


/*     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND */
/*     D INVERSE TO B. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --b;

    /* Function Body */
    k = *n;
L10:
    if (k == 0) {
	goto L80;
    }
    if (kpvt[k] < 0) {
	goto L40;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L30;
    }
    kp = kpvt[k];
    if (kp == k) {
	goto L20;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L20:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    caxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
L30:

/*           APPLY D INVERSE. */

    i__1 = k;
    z_div(&z__1, &b[k], &a[k + k * a_dim1]);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    --k;
    goto L70;
L40:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 2) {
	goto L60;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k - 1) {
	goto L50;
    }

/*                 INTERCHANGE. */

    i__1 = k - 1;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k - 1;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L50:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 2;
    caxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    caxpy_(&i__1, &b[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
L60:

/*           APPLY D INVERSE. */

    z_div(&z__1, &a[k + k * a_dim1], &a[k - 1 + k * a_dim1]);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &a[k - 1 + k * a_dim1]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    z_div(&z__1, &b[k], &a[k - 1 + k * a_dim1]);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &b[k - 1], &a[k - 1 + k * a_dim1]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    k += -2;
L70:
    goto L10;
L80:

/*     LOOP FORWARD APPLYING THE TRANSFORMATIONS. */

    k = 1;
L90:
    if (k > *n) {
	goto L160;
    }
    if (kpvt[k] < 0) {
	goto L120;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L110;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    kp = kpvt[k];
    if (kp == k) {
	goto L100;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L100:
L110:
    ++k;
    goto L150;
L120:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 1) {
	goto L140;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    i__1 = k + 1;
    i__2 = k + 1;
    i__3 = k - 1;
    cdotu_(&z__2, &i__3, &a[(k + 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L130;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L130:
L140:
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* csisl_ */

/* Subroutine */ int csidi_(doublecomplex *a, integer *lda, integer *n, 
	integer *kpvt, doublecomplex *det, doublecomplex *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublecomplex d__;
    static integer j, k;
    static doublecomplex t, ak;
    static integer jb, ks, km1;
    static doublereal ten;
    static doublecomplex akp1, temp, akkp1;
    static logical nodet;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern /* Double Complex */ VOID cdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer kstep;
    static logical noinv;


/*     CSIDI COMPUTES THE DETERMINANT AND INVERSE */
/*     OF A COMPLEX SYMMETRIC MATRIX USING THE FACTORS FROM CSIFA. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA,N) */
/*                THE OUTPUT FROM CSIFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY A. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX A. */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CSIFA. */

/*        WORK    COMPLEX(N) */
/*                WORK VECTOR.  CONTENTS DESTROYED. */

/*        JOB     INTEGER */
/*                JOB HAS THE DECIMAL EXPANSION  AB  WHERE */
/*                   IF  B .NE. 0, THE INVERSE IS COMPUTED, */
/*                   IF  A .NE. 0, THE DETERMINANT IS COMPUTED, */

/*                FOR EXAMPLE, JOB = 11  GIVES BOTH. */

/*     ON RETURN */

/*        VARIABLES NOT REQUESTED BY JOB ARE NOT USED. */

/*        A      CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF */
/*               THE ORIGINAL MATRIX.  THE STRICT LOWER TRIANGLE */
/*               IS NEVER REFERENCED. */

/*        DET    COMPLEX(2) */
/*               DETERMINANT OF ORIGINAL MATRIX. */
/*               DETERMINANT = DET(1) * 10.0**DET(2) */
/*               WITH 1.0 .LE. ABS(DET(1)) .LT. 10.0 */
/*               OR DET(1) = 0.0. */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO MAY OCCUR IF THE INVERSE IS REQUESTED */
/*        AND  CSICO  HAS SET RCOND .EQ. 0.0 */
/*        OR  CSIFA  HAS SET  INFO .NE. 0 . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CCOPY,CDOTU,CSWAP */
/*     FORTRAN ABS,DCMPLX,IABS,MOD,DBLE */

/*     INTERNAL VARIABLES. */



    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --det;
    --work;

    /* Function Body */
    noinv = *job % 10 == 0;
    nodet = *job % 100 / 10 == 0;

    if (nodet) {
	goto L100;
    }
    det[1].r = 1., det[1].i = 0.;
    det[2].r = 0., det[2].i = 0.;
    ten = 10.;
    t.r = 0., t.i = 0.;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * a_dim1;
	d__.r = a[i__2].r, d__.i = a[i__2].i;

/*           CHECK IF 1 BY 1 */

	if (kpvt[k] > 0) {
	    goto L30;
	}

/*              2 BY 2 BLOCK */
/*              USE DET (D  T)  =  (D/T * C - T) * T */
/*                      (T  C) */
/*              TO AVOID UNDERFLOW/OVERFLOW TROUBLES. */
/*              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG. */

	if ((d__1 = t.r, abs(d__1)) + (d__2 = d_imag(&t), abs(d__2)) != 0.) {
	    goto L10;
	}
	i__2 = k + (k + 1) * a_dim1;
	t.r = a[i__2].r, t.i = a[i__2].i;
	z_div(&z__3, &d__, &t);
	i__2 = k + 1 + (k + 1) * a_dim1;
	z__2.r = z__3.r * a[i__2].r - z__3.i * a[i__2].i, z__2.i = z__3.r * a[
		i__2].i + z__3.i * a[i__2].r;
	z__1.r = z__2.r - t.r, z__1.i = z__2.i - t.i;
	d__.r = z__1.r, d__.i = z__1.i;
	goto L20;
L10:
	d__.r = t.r, d__.i = t.i;
	t.r = 0., t.i = 0.;
L20:
L30:

	z__1.r = d__.r * det[1].r - d__.i * det[1].i, z__1.i = d__.r * det[1]
		.i + d__.i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 == 0.) {
	    goto L80;
	}
L40:
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 >= 1.) {
	    goto L50;
	}
	z__2.r = ten, z__2.i = 0.;
	z__1.r = z__2.r * det[1].r - z__2.i * det[1].i, z__1.i = z__2.r * det[
		1].i + z__2.i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r - 1., z__1.i = det[2].i - 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L40;
L50:
L60:
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 < ten) {
	    goto L70;
	}
	z__2.r = ten, z__2.i = 0.;
	z_div(&z__1, &det[1], &z__2);
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r + 1., z__1.i = det[2].i + 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L60;
L70:
L80:
/* L90: */
	;
    }
L100:

/*     COMPUTE INVERSE(A) */

    if (noinv) {
	goto L230;
    }
    k = 1;
L110:
    if (k > *n) {
	goto L220;
    }
    km1 = k - 1;
    if (kpvt[k] < 0) {
	goto L140;
    }

/*              1 BY 1 */

    i__1 = k + k * a_dim1;
    z_div(&z__1, &c_b1092, &a[k + k * a_dim1]);
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L130;
    }
    ccopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + k * a_dim1;
	cdotu_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L120: */
    }
    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    cdotu_(&z__2, &km1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
L130:
    kstep = 1;
    goto L180;
L140:

/*              2 BY 2 */

    i__1 = k + (k + 1) * a_dim1;
    t.r = a[i__1].r, t.i = a[i__1].i;
    z_div(&z__1, &a[k + k * a_dim1], &t);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &t);
    akp1.r = z__1.r, akp1.i = z__1.i;
    z_div(&z__1, &a[k + (k + 1) * a_dim1], &t);
    akkp1.r = z__1.r, akkp1.i = z__1.i;
    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + ak.i * 
	    akp1.r;
    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i * 
	    z__2.r;
    d__.r = z__1.r, d__.i = z__1.i;
    i__1 = k + k * a_dim1;
    z_div(&z__1, &akp1, &d__);
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    i__1 = k + 1 + (k + 1) * a_dim1;
    z_div(&z__1, &ak, &d__);
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    i__1 = k + (k + 1) * a_dim1;
    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
    z_div(&z__1, &z__2, &d__);
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L170;
    }
    ccopy_(&km1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + (k + 1) * a_dim1;
	cdotu_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[(k + 1) * 
		a_dim1 + 1], &c__1);
/* L150: */
    }
    i__1 = k + 1 + (k + 1) * a_dim1;
    i__2 = k + 1 + (k + 1) * a_dim1;
    cdotu_(&z__2, &km1, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    i__1 = k + (k + 1) * a_dim1;
    i__2 = k + (k + 1) * a_dim1;
    cdotu_(&z__2, &km1, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1], &
	    c__1);
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    ccopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + k * a_dim1;
	cdotu_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L160: */
    }
    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    cdotu_(&z__2, &km1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
L170:
    kstep = 2;
L180:

/*           SWAP */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L210;
    }
    cswap_(&ks, &a[ks * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	i__2 = j + k * a_dim1;
	temp.r = a[i__2].r, temp.i = a[i__2].i;
	i__2 = j + k * a_dim1;
	i__3 = ks + j * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = ks + j * a_dim1;
	a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L190: */
    }
    if (kstep == 1) {
	goto L200;
    }
    i__1 = ks + (k + 1) * a_dim1;
    temp.r = a[i__1].r, temp.i = a[i__1].i;
    i__1 = ks + (k + 1) * a_dim1;
    i__2 = k + (k + 1) * a_dim1;
    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
    i__1 = k + (k + 1) * a_dim1;
    a[i__1].r = temp.r, a[i__1].i = temp.i;
L200:
L210:
    k += kstep;
    goto L110;
L220:
L230:
    return 0;
} /* csidi_ */

/* Subroutine */ int cspco_(doublecomplex *ap, integer *n, integer *kpvt, 
	doublereal *rcond, doublecomplex *z__)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static doublecomplex t;
    static integer j1;
    static doublecomplex ak, bk, ek;
    static integer ij, ik, kk, kp, ks, jm1, kps;
    static doublecomplex akm1, bkm1;
    static integer ikm1, km1k, ikp1, info;
    extern /* Subroutine */ int cspfa_(doublecomplex *, integer *, integer *, 
	    integer *);
    static integer km1km1;
    static doublecomplex denom;
    static doublereal anorm;
    extern /* Double Complex */ VOID cdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);


/*     CSPCO FACTORS A COMPLEX SYMMETRIC MATRIX STORED IN PACKED */
/*     FORM BY ELIMINATION WITH SYMMETRIC PIVOTING AND ESTIMATES */
/*     THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, CSPFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW CSPCO BY CSPSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CSPCO BY CSPSL. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW CSPCO BY CSPDI. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CSPCO BY CSPDI. */

/*     ON ENTRY */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE */
/*                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY */
/*                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 . */
/*                SEE COMMENTS BELOW FOR DETAILS. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT STORED IN PACKED FORM. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE */
/*                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        RCOND   REAL */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       COMPLEX(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     PACKED STORAGE */

/*          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER */
/*          TRIANGLE OF A SYMMETRIC MATRIX. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K) = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK CSPFA */
/*     BLAS CAXPY,CDOTU,CSSCAL,SCASUM */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,IABS,DBLE */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A USING ONLY UPPER HALF */

    /* Parameter adjustments */
    --z__;
    --kpvt;
    --ap;

    /* Function Body */
    j1 = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	d__1 = scasum_(&j, &ap[j1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	ij = j1;
	j1 += j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = ij;
	    d__3 = z__[i__4].r + ((d__1 = ap[i__5].r, abs(d__1)) + (d__2 = 
		    d_imag(&ap[ij]), abs(d__2)));
	    z__1.r = d__3, z__1.i = 0.;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    ++ij;
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = j;
	d__1 = anorm, d__2 = z__[i__2].r;
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    cspfa_(&ap[1], n, &kpvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE U*D*W = E */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L50: */
    }
    k = *n;
    ik = *n * (*n - 1) / 2;
L60:
    if (k == 0) {
	goto L120;
    }
    kk = ik + k;
    ikm1 = ik - (k - 1);
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L70;
    }
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L70:
    i__1 = k;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    != 0.) {
	d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	i__2 = k;
	i__3 = k;
	d__8 = (d__5 = z__[i__3].r, abs(d__5)) + (d__6 = d_imag(&z__[k]), abs(
		d__6));
	z__2.r = z__[i__2].r / d__8, z__2.i = z__[i__2].i / d__8;
	z__1.r = d__7 * z__2.r, z__1.i = d__7 * z__2.i;
	ek.r = z__1.r, ek.i = z__1.i;
    }
    i__1 = k;
    i__2 = k;
    z__1.r = z__[i__2].r + ek.r, z__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &ap[ik + 1], &c__1, &z__[1], &c__1);
    if (ks == 1) {
	goto L80;
    }
    i__1 = k - 1;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k - 1]), abs(
	    d__2)) != 0.) {
	d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	i__2 = k - 1;
	i__3 = k - 1;
	d__8 = (d__5 = z__[i__3].r, abs(d__5)) + (d__6 = d_imag(&z__[k - 1]), 
		abs(d__6));
	z__2.r = z__[i__2].r / d__8, z__2.i = z__[i__2].i / d__8;
	z__1.r = d__7 * z__2.r, z__1.i = d__7 * z__2.i;
	ek.r = z__1.r, ek.i = z__1.i;
    }
    i__1 = k - 1;
    i__2 = k - 1;
    z__1.r = z__[i__2].r + ek.r, z__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k - 1], &ap[ikm1 + 1], &c__1, &z__[1], &c__1);
L80:
    if (ks == 2) {
	goto L100;
    }
    i__1 = k;
    i__2 = kk;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    <= (d__3 = ap[i__2].r, abs(d__3)) + (d__4 = d_imag(&ap[kk]), abs(
	    d__4))) {
	goto L90;
    }
    i__1 = kk;
    i__2 = k;
    s = ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)))
	     / ((d__3 = z__[i__2].r, abs(d__3)) + (d__4 = d_imag(&z__[k]), 
	    abs(d__4)));
    csscal_(n, &s, &z__[1], &c__1);
    z__2.r = s, z__2.i = 0.;
    z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + z__2.i * 
	    ek.r;
    ek.r = z__1.r, ek.i = z__1.i;
L90:
    i__1 = kk;
    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)) 
	    != 0.) {
	i__2 = k;
	z_div(&z__1, &z__[k], &ap[kk]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
    }
    i__1 = kk;
    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)) 
	    == 0.) {
	i__2 = k;
	z__[i__2].r = 1., z__[i__2].i = 0.;
    }
    goto L110;
L100:
    km1k = ik + k - 1;
    km1km1 = ikm1 + k - 1;
    z_div(&z__1, &ap[kk], &ap[km1k]);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &ap[km1km1], &ap[km1k]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    z_div(&z__1, &z__[k], &ap[km1k]);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &z__[k - 1], &ap[km1k]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
L110:
    k -= ks;
    ik -= k;
    if (ks == 2) {
	ik -= k + 1;
    }
    goto L60;
L120:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*     SOLVE TRANS(U)*Y = W */

    k = 1;
    ik = 0;
L130:
    if (k > *n) {
	goto L160;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L150;
    }
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&z__2, &i__3, &ap[ik + 1], &c__1, &z__[1], &c__1);
    z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    ikp1 = ik + k;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotu_(&z__2, &i__3, &ap[ikp1 + 1], &c__1, &z__[1], &c__1);
	z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
	z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L140;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L140:
L150:
    ik += k;
    if (ks == 2) {
	ik += k + 1;
    }
    k += ks;
    goto L130;
L160:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE U*D*V = Y */

    k = *n;
    ik = *n * (*n - 1) / 2;
L170:
    if (k == 0) {
	goto L230;
    }
    kk = ik + k;
    ikm1 = ik - (k - 1);
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == ks) {
	goto L190;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L180;
    }
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L180:
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &ap[ik + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - ks;
	caxpy_(&i__1, &z__[k - 1], &ap[ikm1 + 1], &c__1, &z__[1], &c__1);
    }
L190:
    if (ks == 2) {
	goto L210;
    }
    i__1 = k;
    i__2 = kk;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    <= (d__3 = ap[i__2].r, abs(d__3)) + (d__4 = d_imag(&ap[kk]), abs(
	    d__4))) {
	goto L200;
    }
    i__1 = kk;
    i__2 = k;
    s = ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)))
	     / ((d__3 = z__[i__2].r, abs(d__3)) + (d__4 = d_imag(&z__[k]), 
	    abs(d__4)));
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;
L200:
    i__1 = kk;
    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)) 
	    != 0.) {
	i__2 = k;
	z_div(&z__1, &z__[k], &ap[kk]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
    }
    i__1 = kk;
    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)) 
	    == 0.) {
	i__2 = k;
	z__[i__2].r = 1., z__[i__2].i = 0.;
    }
    goto L220;
L210:
    km1k = ik + k - 1;
    km1km1 = ikm1 + k - 1;
    z_div(&z__1, &ap[kk], &ap[km1k]);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &ap[km1km1], &ap[km1k]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    z_div(&z__1, &z__[k], &ap[km1k]);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &z__[k - 1], &ap[km1k]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
L220:
    k -= ks;
    ik -= k;
    if (ks == 2) {
	ik -= k + 1;
    }
    goto L170;
L230:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE TRANS(U)*Z = V */

    k = 1;
    ik = 0;
L240:
    if (k > *n) {
	goto L270;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L260;
    }
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&z__2, &i__3, &ap[ik + 1], &c__1, &z__[1], &c__1);
    z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    ikp1 = ik + k;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotu_(&z__2, &i__3, &ap[ikp1 + 1], &c__1, &z__[1], &c__1);
	z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
	z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L250;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L250:
L260:
    ik += k;
    if (ks == 2) {
	ik += k + 1;
    }
    k += ks;
    goto L240;
L270:
/*     MAKE ZNORM = 1.0 */
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* cspco_ */

/* Subroutine */ int cspfa_(doublecomplex *ap, integer *n, integer *kpvt, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublecomplex t, ak, bk;
    static integer ij, ik, jj, im, jk, kk, km1, km2, ijj, imj, imk;
    static doublecomplex akm1, bkm1;
    static integer ikm1, jkm1, km1k, imim, jmim, imax, jmax;
    static doublecomplex mulk;
    static logical swap;
    static doublereal alpha;
    static integer km1km1;
    static doublecomplex denom;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer kstep, imaxp1;
    static doublecomplex mulkm1;
    static doublereal absakk;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    static doublereal colmax, rowmax;


/*     CSPFA FACTORS A COMPLEX SYMMETRIC MATRIX STORED IN */
/*     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING. */

/*     TO SOLVE  A*X = B , FOLLOW CSPFA BY CSPSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CSPFA BY CSPSL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CSPFA BY CSPDI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW CSPFA BY CSPDI. */

/*     ON ENTRY */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                THE PACKED FORM OF A SYMMETRIC MATRIX  A .  THE */
/*                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY */
/*                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 . */
/*                SEE COMMENTS BELOW FOR DETAILS. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT STORED IN PACKED FORM. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*TRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , TRANS(U) IS THE */
/*                TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS */
/*                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE, */
/*                     BUT IT DOES INDICATE THAT CSPSL OR CSPDI MAY */
/*                     DIVIDE BY ZERO IF CALLED. */

/*     PACKED STORAGE */

/*          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER */
/*          TRIANGLE OF A SYMMETRIC MATRIX. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K)  = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CSWAP,ICAMAX */
/*     FORTRAN ABS,DIMAG,MAX,DBLE,SQRT */

/*     INTERNAL VARIABLES */



/*     INITIALIZE */

/*     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE. */
    /* Parameter adjustments */
    --kpvt;
    --ap;

    /* Function Body */
    alpha = (sqrt(17.) + 1.) / 8.;

    *info = 0;

/*     MAIN LOOP ON K, WHICH GOES FROM N TO 1. */

    k = *n;
    ik = *n * (*n - 1) / 2;
L10:

/*        LEAVE THE LOOP IF K=0 OR K=1. */

/*     ...EXIT */
    if (k == 0) {
	goto L200;
    }
    if (k > 1) {
	goto L20;
    }
    kpvt[1] = 1;
    if ((d__1 = ap[1].r, abs(d__1)) + (d__2 = d_imag(&ap[1]), abs(d__2)) == 
	    0.) {
	*info = 1;
    }
/*     ......EXIT */
    goto L200;
L20:

/*        THIS SECTION OF CODE DETERMINES THE KIND OF */
/*        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED, */
/*        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND */
/*        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS */
/*        REQUIRED. */

    km1 = k - 1;
    kk = ik + k;
    i__1 = kk;
    absakk = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(
	    d__2));

/*        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*        COLUMN K. */

    i__1 = k - 1;
    imax = icamax_(&i__1, &ap[ik + 1], &c__1);
    imk = ik + imax;
    i__1 = imk;
    colmax = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[imk]), abs(
	    d__2));
    if (absakk < alpha * colmax) {
	goto L30;
    }
    kstep = 1;
    swap = FALSE_;
    goto L90;
L30:

/*           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*           ROW IMAX. */

    rowmax = 0.;
    imaxp1 = imax + 1;
    im = imax * (imax - 1) / 2;
    imj = im + (imax << 1);
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = imj;
	d__3 = rowmax, d__4 = (d__1 = ap[i__2].r, abs(d__1)) + (d__2 = d_imag(
		&ap[imj]), abs(d__2));
	rowmax = max(d__3,d__4);
	imj += j;
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = icamax_(&i__1, &ap[im + 1], &c__1);
    jmim = jmax + im;
/* Computing MAX */
    i__1 = jmim;
    d__3 = rowmax, d__4 = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[
	    jmim]), abs(d__2));
    rowmax = max(d__3,d__4);
L50:
    imim = imax + im;
    i__1 = imim;
    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[imim]), abs(d__2))
	     < alpha * rowmax) {
	goto L60;
    }
    kstep = 1;
    swap = TRUE_;
    goto L80;
L60:
    if (absakk < alpha * colmax * (colmax / rowmax)) {
	goto L70;
    }
    kstep = 1;
    swap = FALSE_;
    goto L80;
L70:
    kstep = 2;
    swap = imax != km1;
L80:
L90:
    if (max(absakk,colmax) != 0.) {
	goto L100;
    }

/*           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP. */

    kpvt[k] = k;
    *info = k;
    goto L190;
L100:
    if (kstep == 2) {
	goto L140;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (! swap) {
	goto L120;
    }

/*              PERFORM AN INTERCHANGE. */

    cswap_(&imax, &ap[im + 1], &c__1, &ap[ik + 1], &c__1);
    imj = ik + imax;
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	jk = ik + j;
	i__2 = jk;
	t.r = ap[i__2].r, t.i = ap[i__2].i;
	i__2 = jk;
	i__3 = imj;
	ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
	i__2 = imj;
	ap[i__2].r = t.r, ap[i__2].i = t.i;
	imj -= j - 1;
/* L110: */
    }
L120:

/*           PERFORM THE ELIMINATION. */

    ij = ik - (k - 1);
    i__1 = km1;
    for (jj = 1; jj <= i__1; ++jj) {
	j = k - jj;
	jk = ik + j;
	i__2 = jk;
	z__2.r = -ap[i__2].r, z__2.i = -ap[i__2].i;
	z_div(&z__1, &z__2, &ap[kk]);
	mulk.r = z__1.r, mulk.i = z__1.i;
	t.r = mulk.r, t.i = mulk.i;
	caxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	ijj = ij + j;
	i__2 = jk;
	ap[i__2].r = mulk.r, ap[i__2].i = mulk.i;
	ij -= j - 1;
/* L130: */
    }

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = k;
    if (swap) {
	kpvt[k] = imax;
    }
    goto L190;
L140:

/*           2 X 2 PIVOT BLOCK. */

    km1k = ik + k - 1;
    ikm1 = ik - (k - 1);
    if (! swap) {
	goto L160;
    }

/*              PERFORM AN INTERCHANGE. */

    cswap_(&imax, &ap[im + 1], &c__1, &ap[ikm1 + 1], &c__1);
    imj = ikm1 + imax;
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	jkm1 = ikm1 + j;
	i__2 = jkm1;
	t.r = ap[i__2].r, t.i = ap[i__2].i;
	i__2 = jkm1;
	i__3 = imj;
	ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
	i__2 = imj;
	ap[i__2].r = t.r, ap[i__2].i = t.i;
	imj -= j - 1;
/* L150: */
    }
    i__1 = km1k;
    t.r = ap[i__1].r, t.i = ap[i__1].i;
    i__1 = km1k;
    i__2 = imk;
    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
    i__1 = imk;
    ap[i__1].r = t.r, ap[i__1].i = t.i;
L160:

/*           PERFORM THE ELIMINATION. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    z_div(&z__1, &ap[kk], &ap[km1k]);
    ak.r = z__1.r, ak.i = z__1.i;
    km1km1 = ikm1 + k - 1;
    z_div(&z__1, &ap[km1km1], &ap[km1k]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = 1. - z__2.r, z__1.i = -z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    ij = ik - (k - 1) - (k - 2);
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	jk = ik + j;
	z_div(&z__1, &ap[jk], &ap[km1k]);
	bk.r = z__1.r, bk.i = z__1.i;
	jkm1 = ikm1 + j;
	z_div(&z__1, &ap[jkm1], &ap[km1k]);
	bkm1.r = z__1.r, bkm1.i = z__1.i;
	z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + 
		akm1.i * bk.r;
	z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
	z_div(&z__1, &z__2, &denom);
	mulk.r = z__1.r, mulk.i = z__1.i;
	z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i 
		* bkm1.r;
	z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
	z_div(&z__1, &z__2, &denom);
	mulkm1.r = z__1.r, mulkm1.i = z__1.i;
	t.r = mulk.r, t.i = mulk.i;
	caxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	t.r = mulkm1.r, t.i = mulkm1.i;
	caxpy_(&j, &t, &ap[ikm1 + 1], &c__1, &ap[ij + 1], &c__1);
	i__2 = jk;
	ap[i__2].r = mulk.r, ap[i__2].i = mulk.i;
	i__2 = jkm1;
	ap[i__2].r = mulkm1.r, ap[i__2].i = mulkm1.i;
	ijj = ij + j;
	ij -= j - 1;
/* L170: */
    }
L180:

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = 1 - k;
    if (swap) {
	kpvt[k] = -imax;
    }
    kpvt[k - 1] = kpvt[k];
L190:
    ik -= k - 1;
    if (kstep == 2) {
	ik -= k - 2;
    }
    k -= kstep;
    goto L10;
L200:
    return 0;
} /* cspfa_ */

/* Subroutine */ int cspsl_(doublecomplex *ap, integer *n, integer *kpvt, 
	doublecomplex *b)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex ak, bk;
    static integer ik, kk, kp;
    static doublecomplex akm1, bkm1;
    static integer ikm1, km1k, ikp1;
    static doublecomplex temp;
    static integer km1km1;
    static doublecomplex denom;
    extern /* Double Complex */ VOID cdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CSISL SOLVES THE COMPLEX SYMMETRIC SYSTEM */
/*     A * X = B */
/*     USING THE FACTORS COMPUTED BY CSPFA. */

/*     ON ENTRY */

/*        AP      COMPLEX(N*(N+1)/2) */
/*                THE OUTPUT FROM CSPFA. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CSPFA. */

/*        B       COMPLEX(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO MAY OCCUR IF  CSPCO  HAS SET RCOND .EQ. 0.0 */
/*        OR  CSPFA  HAS SET INFO .NE. 0  . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL CSPFA(AP,N,KPVT,INFO) */
/*           IF (INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CSPSL(AP,N,KPVT,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CDOTU */
/*     FORTRAN IABS */

/*     INTERNAL VARIABLES. */


/*     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND */
/*     D INVERSE TO B. */

    /* Parameter adjustments */
    --b;
    --kpvt;
    --ap;

    /* Function Body */
    k = *n;
    ik = *n * (*n - 1) / 2;
L10:
    if (k == 0) {
	goto L80;
    }
    kk = ik + k;
    if (kpvt[k] < 0) {
	goto L40;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L30;
    }
    kp = kpvt[k];
    if (kp == k) {
	goto L20;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L20:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    caxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
L30:

/*           APPLY D INVERSE. */

    i__1 = k;
    z_div(&z__1, &b[k], &ap[kk]);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    --k;
    ik -= k;
    goto L70;
L40:

/*           2 X 2 PIVOT BLOCK. */

    ikm1 = ik - (k - 1);
    if (k == 2) {
	goto L60;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k - 1) {
	goto L50;
    }

/*                 INTERCHANGE. */

    i__1 = k - 1;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k - 1;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L50:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 2;
    caxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    caxpy_(&i__1, &b[k - 1], &ap[ikm1 + 1], &c__1, &b[1], &c__1);
L60:

/*           APPLY D INVERSE. */

    km1k = ik + k - 1;
    kk = ik + k;
    z_div(&z__1, &ap[kk], &ap[km1k]);
    ak.r = z__1.r, ak.i = z__1.i;
    km1km1 = ikm1 + k - 1;
    z_div(&z__1, &ap[km1km1], &ap[km1k]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    z_div(&z__1, &b[k], &ap[km1k]);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &b[k - 1], &ap[km1k]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    k += -2;
    ik = ik - (k + 1) - k;
L70:
    goto L10;
L80:

/*     LOOP FORWARD APPLYING THE TRANSFORMATIONS. */

    k = 1;
    ik = 0;
L90:
    if (k > *n) {
	goto L160;
    }
    if (kpvt[k] < 0) {
	goto L120;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L110;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&z__2, &i__3, &ap[ik + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    kp = kpvt[k];
    if (kp == k) {
	goto L100;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L100:
L110:
    ik += k;
    ++k;
    goto L150;
L120:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 1) {
	goto L140;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&z__2, &i__3, &ap[ik + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    ikp1 = ik + k;
    i__1 = k + 1;
    i__2 = k + 1;
    i__3 = k - 1;
    cdotu_(&z__2, &i__3, &ap[ikp1 + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L130;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L130:
L140:
    ik = ik + k + k + 1;
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* cspsl_ */

/* Subroutine */ int cspdi_(doublecomplex *ap, integer *n, integer *kpvt, 
	doublecomplex *det, doublecomplex *work, integer *job)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublecomplex d__;
    static integer j, k;
    static doublecomplex t, ak;
    static integer jb, ij, ik, jk, kk, ks, km1;
    static doublereal ten;
    static integer iks, ksj;
    static doublecomplex akp1;
    static integer ikp1, jkp1, kkp1;
    static doublecomplex temp, akkp1;
    static integer kskp1;
    static logical nodet;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern /* Double Complex */ VOID cdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer kstep;
    static logical noinv;


/*     CSPDI COMPUTES THE DETERMINANT AND INVERSE */
/*     OF A COMPLEX SYMMETRIC MATRIX USING THE FACTORS FROM CSPFA, */
/*     WHERE THE MATRIX IS STORED IN PACKED FORM. */

/*     ON ENTRY */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                THE OUTPUT FROM CSPFA. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX A. */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CSPFA. */

/*        WORK    COMPLEX(N) */
/*                WORK VECTOR.  CONTENTS IGNORED. */

/*        JOB     INTEGER */
/*                JOB HAS THE DECIMAL EXPANSION  AB  WHERE */
/*                   IF  B .NE. 0, THE INVERSE IS COMPUTED, */
/*                   IF  A .NE. 0, THE DETERMINANT IS COMPUTED, */

/*                FOR EXAMPLE, JOB = 11  GIVES BOTH. */

/*     ON RETURN */

/*        VARIABLES NOT REQUESTED BY JOB ARE NOT USED. */

/*        AP     CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF */
/*               THE ORIGINAL MATRIX, STORED IN PACKED FORM. */
/*               THE COLUMNS OF THE UPPER TRIANGLE ARE STORED */
/*               SEQUENTIALLY IN A ONE-DIMENSIONAL ARRAY. */

/*        DET    COMPLEX(2) */
/*               DETERMINANT OF ORIGINAL MATRIX. */
/*               DETERMINANT = DET(1) * 10.0**DET(2) */
/*               WITH 1.0 .LE. ABS(DET(1)) .LT. 10.0 */
/*               OR DET(1) = 0.0. */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INVERSE IS REQUESTED */
/*        AND  CSPCO  HAS SET RCOND .EQ. 0.0 */
/*        OR  CSPFA  HAS SET  INFO .NE. 0 . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CCOPY,CDOTU,CSWAP */
/*     FORTRAN ABS,DCMPLX,IABS,MOD,DBLE */

/*     INTERNAL VARIABLES. */



    /* Parameter adjustments */
    --work;
    --det;
    --kpvt;
    --ap;

    /* Function Body */
    noinv = *job % 10 == 0;
    nodet = *job % 100 / 10 == 0;

    if (nodet) {
	goto L110;
    }
    det[1].r = 1., det[1].i = 0.;
    det[2].r = 0., det[2].i = 0.;
    ten = 10.;
    t.r = 0., t.i = 0.;
    ik = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk = ik + k;
	i__2 = kk;
	d__.r = ap[i__2].r, d__.i = ap[i__2].i;

/*           CHECK IF 1 BY 1 */

	if (kpvt[k] > 0) {
	    goto L30;
	}

/*              2 BY 2 BLOCK */
/*              USE DET (D  T)  =  (D/T * C - T) * T */
/*                      (T  C) */
/*              TO AVOID UNDERFLOW/OVERFLOW TROUBLES. */
/*              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG. */

	if ((d__1 = t.r, abs(d__1)) + (d__2 = d_imag(&t), abs(d__2)) != 0.) {
	    goto L10;
	}
	ikp1 = ik + k;
	kkp1 = ikp1 + k;
	i__2 = kkp1;
	t.r = ap[i__2].r, t.i = ap[i__2].i;
	z_div(&z__3, &d__, &t);
	i__2 = kkp1 + 1;
	z__2.r = z__3.r * ap[i__2].r - z__3.i * ap[i__2].i, z__2.i = z__3.r * 
		ap[i__2].i + z__3.i * ap[i__2].r;
	z__1.r = z__2.r - t.r, z__1.i = z__2.i - t.i;
	d__.r = z__1.r, d__.i = z__1.i;
	goto L20;
L10:
	d__.r = t.r, d__.i = t.i;
	t.r = 0., t.i = 0.;
L20:
L30:

	if (nodet) {
	    goto L90;
	}
	z__1.r = d__.r * det[1].r - d__.i * det[1].i, z__1.i = d__.r * det[1]
		.i + d__.i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 == 0.) {
	    goto L80;
	}
L40:
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 >= 1.) {
	    goto L50;
	}
	z__2.r = ten, z__2.i = 0.;
	z__1.r = z__2.r * det[1].r - z__2.i * det[1].i, z__1.i = z__2.r * det[
		1].i + z__2.i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r - 1., z__1.i = det[2].i - 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L40;
L50:
L60:
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 < ten) {
	    goto L70;
	}
	z__2.r = ten, z__2.i = 0.;
	z_div(&z__1, &det[1], &z__2);
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r + 1., z__1.i = det[2].i + 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L60;
L70:
L80:
L90:
	ik += k;
/* L100: */
    }
L110:

/*     COMPUTE INVERSE(A) */

    if (noinv) {
	goto L240;
    }
    k = 1;
    ik = 0;
L120:
    if (k > *n) {
	goto L230;
    }
    km1 = k - 1;
    kk = ik + k;
    ikp1 = ik + k;
    if (kpvt[k] < 0) {
	goto L150;
    }

/*              1 BY 1 */

    i__1 = kk;
    z_div(&z__1, &c_b1092, &ap[kk]);
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L140;
    }
    ccopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	i__2 = jk;
	cdotu_(&z__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L130: */
    }
    i__1 = kk;
    i__2 = kk;
    cdotu_(&z__2, &km1, &work[1], &c__1, &ap[ik + 1], &c__1);
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
L140:
    kstep = 1;
    goto L190;
L150:

/*              2 BY 2 */

    kkp1 = ikp1 + k;
    i__1 = kkp1;
    t.r = ap[i__1].r, t.i = ap[i__1].i;
    z_div(&z__1, &ap[kk], &t);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &ap[kkp1 + 1], &t);
    akp1.r = z__1.r, akp1.i = z__1.i;
    z_div(&z__1, &ap[kkp1], &t);
    akkp1.r = z__1.r, akkp1.i = z__1.i;
    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + ak.i * 
	    akp1.r;
    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i * 
	    z__2.r;
    d__.r = z__1.r, d__.i = z__1.i;
    i__1 = kk;
    z_div(&z__1, &akp1, &d__);
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    i__1 = kkp1 + 1;
    z_div(&z__1, &ak, &d__);
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    i__1 = kkp1;
    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
    z_div(&z__1, &z__2, &d__);
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L180;
    }
    ccopy_(&km1, &ap[ikp1 + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jkp1 = ikp1 + j;
	i__2 = jkp1;
	cdotu_(&z__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ikp1 + 1], &c__1);
	ij += j;
/* L160: */
    }
    i__1 = kkp1 + 1;
    i__2 = kkp1 + 1;
    cdotu_(&z__2, &km1, &work[1], &c__1, &ap[ikp1 + 1], &c__1);
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    i__1 = kkp1;
    i__2 = kkp1;
    cdotu_(&z__2, &km1, &ap[ik + 1], &c__1, &ap[ikp1 + 1], &c__1);
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    ccopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	i__2 = jk;
	cdotu_(&z__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L170: */
    }
    i__1 = kk;
    i__2 = kk;
    cdotu_(&z__2, &km1, &work[1], &c__1, &ap[ik + 1], &c__1);
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
L180:
    kstep = 2;
L190:

/*           SWAP */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L220;
    }
    iks = ks * (ks - 1) / 2;
    cswap_(&ks, &ap[iks + 1], &c__1, &ap[ik + 1], &c__1);
    ksj = ik + ks;
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	jk = ik + j;
	i__2 = jk;
	temp.r = ap[i__2].r, temp.i = ap[i__2].i;
	i__2 = jk;
	i__3 = ksj;
	ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
	i__2 = ksj;
	ap[i__2].r = temp.r, ap[i__2].i = temp.i;
	ksj -= j - 1;
/* L200: */
    }
    if (kstep == 1) {
	goto L210;
    }
    kskp1 = ikp1 + ks;
    i__1 = kskp1;
    temp.r = ap[i__1].r, temp.i = ap[i__1].i;
    i__1 = kskp1;
    i__2 = kkp1;
    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
    i__1 = kkp1;
    ap[i__1].r = temp.r, ap[i__1].i = temp.i;
L210:
L220:
    ik += k;
    if (kstep == 2) {
	ik = ik + k + 1;
    }
    k += kstep;
    goto L120;
L230:
L240:
    return 0;
} /* cspdi_ */

/* Subroutine */ int chico_(doublecomplex *a, integer *lda, integer *n, 
	integer *kpvt, doublereal *rcond, doublecomplex *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static doublecomplex t, ak, bk, ek;
    static integer kp, ks, jm1, kps;
    static doublecomplex akm1, bkm1;
    static integer info;
    extern /* Subroutine */ int chifa_(doublecomplex *, integer *, integer *, 
	    integer *, integer *);
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublecomplex denom;
    static doublereal anorm;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);


/*     CHICO FACTORS A COMPLEX HERMITIAN MATRIX BY ELIMINATION WITH */
/*     SYMMETRIC PIVOTING AND ESTIMATES THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, CHIFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW CHICO BY CHISL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CHICO BY CHISL. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW CHICO BY CHIDI. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CHICO BY CHIDI. */
/*     TO COMPUTE  INERTIA(A), FOLLOW CHICO BY CHIDI. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA, N) */
/*                THE HERMITIAN MATRIX TO BE FACTORED. */
/*                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     OUTPUT */

/*        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*CTRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , CTRANS(U) IS THE */
/*                CONJUGATE TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        RCOND   REAL */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       COMPLEX(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK CHIFA */
/*     BLAS CAXPY,CDOTC,CSSCAL,SCASUM */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,CONJG,IABS,DBLE */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A USING ONLY UPPER HALF */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	d__1 = scasum_(&j, &a[j * a_dim1 + 1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__ + j * a_dim1;
	    d__3 = z__[i__4].r + ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = 
		    d_imag(&a[i__ + j * a_dim1]), abs(d__2)));
	    z__1.r = d__3, z__1.i = 0.;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = j;
	d__1 = anorm, d__2 = z__[i__2].r;
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    chifa_(&a[a_offset], lda, n, &kpvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE U*D*W = E */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L50: */
    }
    k = *n;
L60:
    if (k == 0) {
	goto L120;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L70;
    }
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L70:
    i__1 = k;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    != 0.) {
	d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	i__2 = k;
	i__3 = k;
	d__8 = (d__5 = z__[i__3].r, abs(d__5)) + (d__6 = d_imag(&z__[k]), abs(
		d__6));
	z__2.r = z__[i__2].r / d__8, z__2.i = z__[i__2].i / d__8;
	z__1.r = d__7 * z__2.r, z__1.i = d__7 * z__2.i;
	ek.r = z__1.r, ek.i = z__1.i;
    }
    i__1 = k;
    i__2 = k;
    z__1.r = z__[i__2].r + ek.r, z__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 1) {
	goto L80;
    }
    i__1 = k - 1;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k - 1]), abs(
	    d__2)) != 0.) {
	d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	i__2 = k - 1;
	i__3 = k - 1;
	d__8 = (d__5 = z__[i__3].r, abs(d__5)) + (d__6 = d_imag(&z__[k - 1]), 
		abs(d__6));
	z__2.r = z__[i__2].r / d__8, z__2.i = z__[i__2].i / d__8;
	z__1.r = d__7 * z__2.r, z__1.i = d__7 * z__2.i;
	ek.r = z__1.r, ek.i = z__1.i;
    }
    i__1 = k - 1;
    i__2 = k - 1;
    z__1.r = z__[i__2].r + ek.r, z__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
	    c__1);
L80:
    if (ks == 2) {
	goto L100;
    }
    i__1 = k;
    i__2 = k + k * a_dim1;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    <= (d__3 = a[i__2].r, abs(d__3)) + (d__4 = d_imag(&a[k + k * 
	    a_dim1]), abs(d__4))) {
	goto L90;
    }
    i__1 = k + k * a_dim1;
    i__2 = k;
    s = ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2))) / ((d__3 = z__[i__2].r, abs(d__3)) + (d__4 = d_imag(&
	    z__[k]), abs(d__4)));
    csscal_(n, &s, &z__[1], &c__1);
    z__2.r = s, z__2.i = 0.;
    z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + z__2.i * 
	    ek.r;
    ek.r = z__1.r, ek.i = z__1.i;
L90:
    i__1 = k + k * a_dim1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2)) != 0.) {
	i__2 = k;
	z_div(&z__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
    }
    i__1 = k + k * a_dim1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2)) == 0.) {
	i__2 = k;
	z__[i__2].r = 1., z__[i__2].i = 0.;
    }
    goto L110;
L100:
    d_cnjg(&z__2, &a[k - 1 + k * a_dim1]);
    z_div(&z__1, &a[k + k * a_dim1], &z__2);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &a[k - 1 + k * a_dim1]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    d_cnjg(&z__2, &a[k - 1 + k * a_dim1]);
    z_div(&z__1, &z__[k], &z__2);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &z__[k - 1], &a[k - 1 + k * a_dim1]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
L110:
    k -= ks;
    goto L60;
L120:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*     SOLVE CTRANS(U)*Y = W */

    k = 1;
L130:
    if (k > *n) {
	goto L160;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L150;
    }
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotc_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotc_(&z__2, &i__3, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &c__1);
	z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
	z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L140;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L140:
L150:
    k += ks;
    goto L130;
L160:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE U*D*V = Y */

    k = *n;
L170:
    if (k == 0) {
	goto L230;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == ks) {
	goto L190;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L180;
    }
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L180:
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - ks;
	caxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
L190:
    if (ks == 2) {
	goto L210;
    }
    i__1 = k;
    i__2 = k + k * a_dim1;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    <= (d__3 = a[i__2].r, abs(d__3)) + (d__4 = d_imag(&a[k + k * 
	    a_dim1]), abs(d__4))) {
	goto L200;
    }
    i__1 = k + k * a_dim1;
    i__2 = k;
    s = ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2))) / ((d__3 = z__[i__2].r, abs(d__3)) + (d__4 = d_imag(&
	    z__[k]), abs(d__4)));
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;
L200:
    i__1 = k + k * a_dim1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2)) != 0.) {
	i__2 = k;
	z_div(&z__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
    }
    i__1 = k + k * a_dim1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]), 
	    abs(d__2)) == 0.) {
	i__2 = k;
	z__[i__2].r = 1., z__[i__2].i = 0.;
    }
    goto L220;
L210:
    d_cnjg(&z__2, &a[k - 1 + k * a_dim1]);
    z_div(&z__1, &a[k + k * a_dim1], &z__2);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &a[k - 1 + k * a_dim1]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    d_cnjg(&z__2, &a[k - 1 + k * a_dim1]);
    z_div(&z__1, &z__[k], &z__2);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &z__[k - 1], &a[k - 1 + k * a_dim1]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
L220:
    k -= ks;
    goto L170;
L230:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE CTRANS(U)*Z = V */

    k = 1;
L240:
    if (k > *n) {
	goto L270;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L260;
    }
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotc_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotc_(&z__2, &i__3, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &c__1);
	z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
	z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L250;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L250:
L260:
    k += ks;
    goto L240;
L270:
/*     MAKE ZNORM = 1.0 */
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* chico_ */

/* Subroutine */ int chifa_(doublecomplex *a, integer *lda, integer *n, 
	integer *kpvt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublecomplex t, ak, bk;
    static integer jj, km1, km2;
    static doublecomplex akm1, bkm1;
    static integer imax, jmax;
    static doublecomplex mulk;
    static logical swap;
    static doublereal alpha;
    static doublecomplex denom;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer kstep, imaxp1;
    static doublecomplex mulkm1;
    static doublereal absakk;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    static doublereal colmax, rowmax;


/*     CHIFA FACTORS A COMPLEX HERMITIAN MATRIX BY ELIMINATION */
/*     WITH SYMMETRIC PIVOTING. */

/*     TO SOLVE  A*X = B , FOLLOW CHIFA BY CHISL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CHIFA BY CHISL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CHIFA BY CHIDI. */
/*     TO COMPUTE  INERTIA(A) , FOLLOW CHIFA BY CHIDI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW CHIFA BY CHIDI. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA,N) */
/*                THE HERMITIAN MATRIX TO BE FACTORED. */
/*                ONLY THE DIAGONAL AND UPPER TRIANGLE ARE USED. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     ON RETURN */

/*        A       A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*CTRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , CTRANS(U) IS THE */
/*                CONJUGATE TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS */
/*                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE, */
/*                     BUT IT DOES INDICATE THAT CHISL OR CHIDI MAY */
/*                     DIVIDE BY ZERO IF CALLED. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CSWAP,ICAMAX */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,CONJG,DBLE,SQRT */

/*     INTERNAL VARIABLES */



/*     INITIALIZE */

/*     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;

    /* Function Body */
    alpha = (sqrt(17.) + 1.) / 8.;

    *info = 0;

/*     MAIN LOOP ON K, WHICH GOES FROM N TO 1. */

    k = *n;
L10:

/*        LEAVE THE LOOP IF K=0 OR K=1. */

/*     ...EXIT */
    if (k == 0) {
	goto L200;
    }
    if (k > 1) {
	goto L20;
    }
    kpvt[1] = 1;
    i__1 = a_dim1 + 1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[a_dim1 + 1]), abs(
	    d__2)) == 0.) {
	*info = 1;
    }
/*     ......EXIT */
    goto L200;
L20:

/*        THIS SECTION OF CODE DETERMINES THE KIND OF */
/*        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED, */
/*        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND */
/*        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS */
/*        REQUIRED. */

    km1 = k - 1;
    i__1 = k + k * a_dim1;
    absakk = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[k + k * a_dim1]
	    ), abs(d__2));

/*        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*        COLUMN K. */

    i__1 = k - 1;
    imax = icamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = imax + k * a_dim1;
    colmax = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + k * 
	    a_dim1]), abs(d__2));
    if (absakk < alpha * colmax) {
	goto L30;
    }
    kstep = 1;
    swap = FALSE_;
    goto L90;
L30:

/*           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*           ROW IMAX. */

    rowmax = 0.;
    imaxp1 = imax + 1;
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = imax + j * a_dim1;
	d__3 = rowmax, d__4 = (d__1 = a[i__2].r, abs(d__1)) + (d__2 = d_imag(&
		a[imax + j * a_dim1]), abs(d__2));
	rowmax = max(d__3,d__4);
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = icamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
    i__1 = jmax + imax * a_dim1;
    d__3 = rowmax, d__4 = (d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[
	    jmax + imax * a_dim1]), abs(d__2));
    rowmax = max(d__3,d__4);
L50:
    i__1 = imax + imax * a_dim1;
    if ((d__1 = a[i__1].r, abs(d__1)) + (d__2 = d_imag(&a[imax + imax * 
	    a_dim1]), abs(d__2)) < alpha * rowmax) {
	goto L60;
    }
    kstep = 1;
    swap = TRUE_;
    goto L80;
L60:
    if (absakk < alpha * colmax * (colmax / rowmax)) {
	goto L70;
    }
    kstep = 1;
    swap = FALSE_;
    goto L80;
L70:
    kstep = 2;
    swap = imax != km1;
L80:
L90:
    if (max(absakk,colmax) != 0.) {
	goto L100;
    }

/*           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP. */

    kpvt[k] = k;
    *info = k;
    goto L190;
L100:
    if (kstep == 2) {
	goto L140;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (! swap) {
	goto L120;
    }

/*              PERFORM AN INTERCHANGE. */

    cswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	d_cnjg(&z__1, &a[j + k * a_dim1]);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = j + k * a_dim1;
	d_cnjg(&z__1, &a[imax + j * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = imax + j * a_dim1;
	a[i__2].r = t.r, a[i__2].i = t.i;
/* L110: */
    }
L120:

/*           PERFORM THE ELIMINATION. */

    i__1 = km1;
    for (jj = 1; jj <= i__1; ++jj) {
	j = k - jj;
	i__2 = j + k * a_dim1;
	z__2.r = -a[i__2].r, z__2.i = -a[i__2].i;
	z_div(&z__1, &z__2, &a[k + k * a_dim1]);
	mulk.r = z__1.r, mulk.i = z__1.i;
	d_cnjg(&z__1, &mulk);
	t.r = z__1.r, t.i = z__1.i;
	caxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	i__2 = j + j * a_dim1;
	i__3 = j + j * a_dim1;
	d__1 = a[i__3].r;
	z__1.r = d__1, z__1.i = 0.;
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j + k * a_dim1;
	a[i__2].r = mulk.r, a[i__2].i = mulk.i;
/* L130: */
    }

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = k;
    if (swap) {
	kpvt[k] = imax;
    }
    goto L190;
L140:

/*           2 X 2 PIVOT BLOCK. */

    if (! swap) {
	goto L160;
    }

/*              PERFORM AN INTERCHANGE. */

    cswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[(k - 1) * a_dim1 + 1], &
	    c__1);
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	d_cnjg(&z__1, &a[j + (k - 1) * a_dim1]);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = j + (k - 1) * a_dim1;
	d_cnjg(&z__1, &a[imax + j * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = imax + j * a_dim1;
	a[i__2].r = t.r, a[i__2].i = t.i;
/* L150: */
    }
    i__1 = k - 1 + k * a_dim1;
    t.r = a[i__1].r, t.i = a[i__1].i;
    i__1 = k - 1 + k * a_dim1;
    i__2 = imax + k * a_dim1;
    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
    i__1 = imax + k * a_dim1;
    a[i__1].r = t.r, a[i__1].i = t.i;
L160:

/*           PERFORM THE ELIMINATION. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    z_div(&z__1, &a[k + k * a_dim1], &a[k - 1 + k * a_dim1]);
    ak.r = z__1.r, ak.i = z__1.i;
    d_cnjg(&z__2, &a[k - 1 + k * a_dim1]);
    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &z__2);
    akm1.r = z__1.r, akm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = 1. - z__2.r, z__1.i = -z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	z_div(&z__1, &a[j + k * a_dim1], &a[k - 1 + k * a_dim1]);
	bk.r = z__1.r, bk.i = z__1.i;
	d_cnjg(&z__2, &a[k - 1 + k * a_dim1]);
	z_div(&z__1, &a[j + (k - 1) * a_dim1], &z__2);
	bkm1.r = z__1.r, bkm1.i = z__1.i;
	z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + 
		akm1.i * bk.r;
	z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
	z_div(&z__1, &z__2, &denom);
	mulk.r = z__1.r, mulk.i = z__1.i;
	z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i 
		* bkm1.r;
	z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
	z_div(&z__1, &z__2, &denom);
	mulkm1.r = z__1.r, mulkm1.i = z__1.i;
	d_cnjg(&z__1, &mulk);
	t.r = z__1.r, t.i = z__1.i;
	caxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	d_cnjg(&z__1, &mulkm1);
	t.r = z__1.r, t.i = z__1.i;
	caxpy_(&j, &t, &a[(k - 1) * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		c__1);
	i__2 = j + k * a_dim1;
	a[i__2].r = mulk.r, a[i__2].i = mulk.i;
	i__2 = j + (k - 1) * a_dim1;
	a[i__2].r = mulkm1.r, a[i__2].i = mulkm1.i;
	i__2 = j + j * a_dim1;
	i__3 = j + j * a_dim1;
	d__1 = a[i__3].r;
	z__1.r = d__1, z__1.i = 0.;
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
/* L170: */
    }
L180:

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = 1 - k;
    if (swap) {
	kpvt[k] = -imax;
    }
    kpvt[k - 1] = kpvt[k];
L190:
    k -= kstep;
    goto L10;
L200:
    return 0;
} /* chifa_ */

/* Subroutine */ int chisl_(doublecomplex *a, integer *lda, integer *n, 
	integer *kpvt, doublecomplex *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex ak, bk;
    static integer kp;
    static doublecomplex akm1, bkm1, temp;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublecomplex denom;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CHISL SOLVES THE COMPLEX HERMITIAN SYSTEM */
/*     A * X = B */
/*     USING THE FACTORS COMPUTED BY CHIFA. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA,N) */
/*                THE OUTPUT FROM CHIFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY  A . */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CHIFA. */

/*        B       COMPLEX(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO MAY OCCUR IF  CHICO  HAS SET RCOND .EQ. 0.0 */
/*        OR  CHIFA  HAS SET INFO .NE. 0  . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL CHIFA(A,LDA,N,KPVT,INFO) */
/*           IF (INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CHISL(A,LDA,N,KPVT,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CDOTC */
/*     FORTRAN CONJG,IABS */

/*     INTERNAL VARIABLES. */


/*     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND */
/*     D INVERSE TO B. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --b;

    /* Function Body */
    k = *n;
L10:
    if (k == 0) {
	goto L80;
    }
    if (kpvt[k] < 0) {
	goto L40;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L30;
    }
    kp = kpvt[k];
    if (kp == k) {
	goto L20;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L20:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    caxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
L30:

/*           APPLY D INVERSE. */

    i__1 = k;
    z_div(&z__1, &b[k], &a[k + k * a_dim1]);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    --k;
    goto L70;
L40:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 2) {
	goto L60;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k - 1) {
	goto L50;
    }

/*                 INTERCHANGE. */

    i__1 = k - 1;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k - 1;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L50:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 2;
    caxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    caxpy_(&i__1, &b[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
L60:

/*           APPLY D INVERSE. */

    d_cnjg(&z__2, &a[k - 1 + k * a_dim1]);
    z_div(&z__1, &a[k + k * a_dim1], &z__2);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &a[k - 1 + (k - 1) * a_dim1], &a[k - 1 + k * a_dim1]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    d_cnjg(&z__2, &a[k - 1 + k * a_dim1]);
    z_div(&z__1, &b[k], &z__2);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &b[k - 1], &a[k - 1 + k * a_dim1]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    k += -2;
L70:
    goto L10;
L80:

/*     LOOP FORWARD APPLYING THE TRANSFORMATIONS. */

    k = 1;
L90:
    if (k > *n) {
	goto L160;
    }
    if (kpvt[k] < 0) {
	goto L120;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L110;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotc_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    kp = kpvt[k];
    if (kp == k) {
	goto L100;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L100:
L110:
    ++k;
    goto L150;
L120:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 1) {
	goto L140;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotc_(&z__2, &i__3, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    i__1 = k + 1;
    i__2 = k + 1;
    i__3 = k - 1;
    cdotc_(&z__2, &i__3, &a[(k + 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L130;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L130:
L140:
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* chisl_ */

/* Subroutine */ int chidi_(doublecomplex *a, integer *lda, integer *n, 
	integer *kpvt, doublereal *det, integer *inert, doublecomplex *work, 
	integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal d__;
    static integer j, k;
    static doublereal t, ak;
    static integer jb, ks, km1;
    static doublereal ten, akp1;
    static doublecomplex temp, akkp1;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical nodet;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), caxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    static integer kstep;
    static logical noert, noinv;


/*     CHIDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE */
/*     OF A COMPLEX HERMITIAN MATRIX USING THE FACTORS FROM CHIFA. */

/*     ON ENTRY */

/*        A       COMPLEX(LDA,N) */
/*                THE OUTPUT FROM CHIFA. */

/*        LDA     INTEGER */
/*                THE LEADING DIMENSION OF THE ARRAY A. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX A. */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CHIFA. */

/*        WORK    COMPLEX(N) */
/*                WORK VECTOR.  CONTENTS DESTROYED. */

/*        JOB     INTEGER */
/*                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE */
/*                   IF  C .NE. 0, THE INVERSE IS COMPUTED, */
/*                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED, */
/*                   IF  A .NE. 0, THE INERTIA IS COMPUTED. */

/*                FOR EXAMPLE, JOB = 111  GIVES ALL THREE. */

/*     ON RETURN */

/*        VARIABLES NOT REQUESTED BY JOB ARE NOT USED. */

/*        A      CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF */
/*               THE ORIGINAL MATRIX.  THE STRICT LOWER TRIANGLE */
/*               IS NEVER REFERENCED. */

/*        DET    REAL(2) */
/*               DETERMINANT OF ORIGINAL MATRIX. */
/*               DETERMINANT = DET(1) * 10.0**DET(2) */
/*               WITH 1.0 .LE. ABS(DET(1)) .LT. 10.0 */
/*               OR DET(1) = 0.0. */

/*        INERT  INTEGER(3) */
/*               THE INERTIA OF THE ORIGINAL MATRIX. */
/*               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES. */
/*               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES. */
/*               INERT(3)  =  NUMBER OF ZERO EIGENVALUES. */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO MAY OCCUR IF THE INVERSE IS REQUESTED */
/*        AND  CHICO  HAS SET RCOND .EQ. 0.0 */
/*        OR  CHIFA  HAS SET  INFO .NE. 0 . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CCOPY,CDOTC,CSWAP */
/*     FORTRAN ABS,CABS,DCMPLX,CONJG,IABS,MOD,DBLE */

/*     INTERNAL VARIABLES. */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --det;
    --inert;
    --work;

    /* Function Body */
    noinv = *job % 10 == 0;
    nodet = *job % 100 / 10 == 0;
    noert = *job % 1000 / 100 == 0;

    if (nodet && noert) {
	goto L140;
    }
    if (noert) {
	goto L10;
    }
    inert[1] = 0;
    inert[2] = 0;
    inert[3] = 0;
L10:
    if (nodet) {
	goto L20;
    }
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
L20:
    t = 0.;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * a_dim1;
	d__ = a[i__2].r;

/*           CHECK IF 1 BY 1 */

	if (kpvt[k] > 0) {
	    goto L50;
	}

/*              2 BY 2 BLOCK */
/*              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = ABS(S) */
/*                      (S  C) */
/*              TO AVOID UNDERFLOW/OVERFLOW TROUBLES. */
/*              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG. */

	if (t != 0.) {
	    goto L30;
	}
	t = z_abs(&a[k + (k + 1) * a_dim1]);
	i__2 = k + 1 + (k + 1) * a_dim1;
	d__ = d__ / t * a[i__2].r - t;
	goto L40;
L30:
	d__ = t;
	t = 0.;
L40:
L50:

	if (noert) {
	    goto L60;
	}
	if (d__ > 0.) {
	    ++inert[1];
	}
	if (d__ < 0.) {
	    ++inert[2];
	}
	if (d__ == 0.) {
	    ++inert[3];
	}
L60:

	if (nodet) {
	    goto L120;
	}
	det[1] = d__ * det[1];
	if (det[1] == 0.) {
	    goto L110;
	}
L70:
	if (abs(det[1]) >= 1.) {
	    goto L80;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L70;
L80:
L90:
	if (abs(det[1]) < ten) {
	    goto L100;
	}
	det[1] /= ten;
	det[2] += 1.;
	goto L90;
L100:
L110:
L120:
/* L130: */
	;
    }
L140:

/*     COMPUTE INVERSE(A) */

    if (noinv) {
	goto L270;
    }
    k = 1;
L150:
    if (k > *n) {
	goto L260;
    }
    km1 = k - 1;
    if (kpvt[k] < 0) {
	goto L180;
    }

/*              1 BY 1 */

    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    d__1 = 1. / a[i__2].r;
    z__1.r = d__1, z__1.i = 0.;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L170;
    }
    ccopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + k * a_dim1;
	cdotc_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L160: */
    }
    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    cdotc_(&z__3, &km1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
    d__1 = z__3.r;
    z__2.r = d__1, z__2.i = 0.;
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
L170:
    kstep = 1;
    goto L220;
L180:

/*              2 BY 2 */

    t = z_abs(&a[k + (k + 1) * a_dim1]);
    i__1 = k + k * a_dim1;
    ak = a[i__1].r / t;
    i__1 = k + 1 + (k + 1) * a_dim1;
    akp1 = a[i__1].r / t;
    i__1 = k + (k + 1) * a_dim1;
    z__1.r = a[i__1].r / t, z__1.i = a[i__1].i / t;
    akkp1.r = z__1.r, akkp1.i = z__1.i;
    d__ = t * (ak * akp1 - 1.);
    i__1 = k + k * a_dim1;
    d__1 = akp1 / d__;
    z__1.r = d__1, z__1.i = 0.;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    i__1 = k + 1 + (k + 1) * a_dim1;
    d__1 = ak / d__;
    z__1.r = d__1, z__1.i = 0.;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    i__1 = k + (k + 1) * a_dim1;
    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
    z__1.r = z__2.r / d__, z__1.i = z__2.i / d__;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L210;
    }
    ccopy_(&km1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + (k + 1) * a_dim1;
	cdotc_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[(k + 1) * 
		a_dim1 + 1], &c__1);
/* L190: */
    }
    i__1 = k + 1 + (k + 1) * a_dim1;
    i__2 = k + 1 + (k + 1) * a_dim1;
    cdotc_(&z__3, &km1, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
    d__1 = z__3.r;
    z__2.r = d__1, z__2.i = 0.;
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    i__1 = k + (k + 1) * a_dim1;
    i__2 = k + (k + 1) * a_dim1;
    cdotc_(&z__2, &km1, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1], &
	    c__1);
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    ccopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + k * a_dim1;
	cdotc_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L200: */
    }
    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    cdotc_(&z__3, &km1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
    d__1 = z__3.r;
    z__2.r = d__1, z__2.i = 0.;
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
L210:
    kstep = 2;
L220:

/*           SWAP */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L250;
    }
    cswap_(&ks, &a[ks * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	d_cnjg(&z__1, &a[j + k * a_dim1]);
	temp.r = z__1.r, temp.i = z__1.i;
	i__2 = j + k * a_dim1;
	d_cnjg(&z__1, &a[ks + j * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = ks + j * a_dim1;
	a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L230: */
    }
    if (kstep == 1) {
	goto L240;
    }
    i__1 = ks + (k + 1) * a_dim1;
    temp.r = a[i__1].r, temp.i = a[i__1].i;
    i__1 = ks + (k + 1) * a_dim1;
    i__2 = k + (k + 1) * a_dim1;
    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
    i__1 = k + (k + 1) * a_dim1;
    a[i__1].r = temp.r, a[i__1].i = temp.i;
L240:
L250:
    k += kstep;
    goto L150;
L260:
L270:
    return 0;
} /* chidi_ */

/* Subroutine */ int chpco_(doublecomplex *ap, integer *n, integer *kpvt, 
	doublereal *rcond, doublecomplex *z__)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal s;
    static doublecomplex t;
    static integer j1;
    static doublecomplex ak, bk, ek;
    static integer ij, ik, kk, kp, ks, jm1, kps;
    static doublecomplex akm1, bkm1;
    static integer ikm1, km1k, ikp1, info;
    extern /* Subroutine */ int chpfa_(doublecomplex *, integer *, integer *, 
	    integer *);
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer km1km1;
    static doublecomplex denom;
    static doublereal anorm;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal ynorm;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);


/*     CHPCO FACTORS A COMPLEX HERMITIAN MATRIX STORED IN PACKED */
/*     FORM BY ELIMINATION WITH SYMMETRIC PIVOTING AND ESTIMATES */
/*     THE CONDITION OF THE MATRIX. */

/*     IF  RCOND  IS NOT NEEDED, CHPFA IS SLIGHTLY FASTER. */
/*     TO SOLVE  A*X = B , FOLLOW CHPCO BY CHPSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CHPCO BY CHPSL. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW CHPCO BY CHPDI. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CHPCO BY CHPDI. */
/*     TO COMPUTE  INERTIA(A), FOLLOW CHPCO BY CHPDI. */

/*     ON ENTRY */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                THE PACKED FORM OF A HERMITIAN MATRIX  A .  THE */
/*                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY */
/*                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 . */
/*                SEE COMMENTS BELOW FOR DETAILS. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     OUTPUT */

/*        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT STORED IN PACKED FORM. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*CTRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , CTRANS(U) IS THE */
/*                CONJUGATE TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        RCOND   REAL */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A . */
/*                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS */
/*                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       COMPLEX(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     PACKED STORAGE */

/*          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER */
/*          TRIANGLE OF A HERMITIAN MATRIX. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K) = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     LINPACK CHPFA */
/*     BLAS CAXPY,CDOTC,CSSCAL,SCASUM */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,CONJG,IABS,DBLE */

/*     INTERNAL VARIABLES */



/*     FIND NORM OF A USING ONLY UPPER HALF */

    /* Parameter adjustments */
    --z__;
    --kpvt;
    --ap;

    /* Function Body */
    j1 = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	d__1 = scasum_(&j, &ap[j1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	ij = j1;
	j1 += j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = ij;
	    d__3 = z__[i__4].r + ((d__1 = ap[i__5].r, abs(d__1)) + (d__2 = 
		    d_imag(&ap[ij]), abs(d__2)));
	    z__1.r = d__3, z__1.i = 0.;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    ++ij;
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = j;
	d__1 = anorm, d__2 = z__[i__2].r;
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     FACTOR */

    chpfa_(&ap[1], n, &kpvt[1], &info);

/*     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF W  WHERE  U*D*W = E . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE U*D*W = E */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L50: */
    }
    k = *n;
    ik = *n * (*n - 1) / 2;
L60:
    if (k == 0) {
	goto L120;
    }
    kk = ik + k;
    ikm1 = ik - (k - 1);
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L70;
    }
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L70:
    i__1 = k;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    != 0.) {
	d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	i__2 = k;
	i__3 = k;
	d__8 = (d__5 = z__[i__3].r, abs(d__5)) + (d__6 = d_imag(&z__[k]), abs(
		d__6));
	z__2.r = z__[i__2].r / d__8, z__2.i = z__[i__2].i / d__8;
	z__1.r = d__7 * z__2.r, z__1.i = d__7 * z__2.i;
	ek.r = z__1.r, ek.i = z__1.i;
    }
    i__1 = k;
    i__2 = k;
    z__1.r = z__[i__2].r + ek.r, z__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &ap[ik + 1], &c__1, &z__[1], &c__1);
    if (ks == 1) {
	goto L80;
    }
    i__1 = k - 1;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k - 1]), abs(
	    d__2)) != 0.) {
	d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	i__2 = k - 1;
	i__3 = k - 1;
	d__8 = (d__5 = z__[i__3].r, abs(d__5)) + (d__6 = d_imag(&z__[k - 1]), 
		abs(d__6));
	z__2.r = z__[i__2].r / d__8, z__2.i = z__[i__2].i / d__8;
	z__1.r = d__7 * z__2.r, z__1.i = d__7 * z__2.i;
	ek.r = z__1.r, ek.i = z__1.i;
    }
    i__1 = k - 1;
    i__2 = k - 1;
    z__1.r = z__[i__2].r + ek.r, z__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k - 1], &ap[ikm1 + 1], &c__1, &z__[1], &c__1);
L80:
    if (ks == 2) {
	goto L100;
    }
    i__1 = k;
    i__2 = kk;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    <= (d__3 = ap[i__2].r, abs(d__3)) + (d__4 = d_imag(&ap[kk]), abs(
	    d__4))) {
	goto L90;
    }
    i__1 = kk;
    i__2 = k;
    s = ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)))
	     / ((d__3 = z__[i__2].r, abs(d__3)) + (d__4 = d_imag(&z__[k]), 
	    abs(d__4)));
    csscal_(n, &s, &z__[1], &c__1);
    z__2.r = s, z__2.i = 0.;
    z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + z__2.i * 
	    ek.r;
    ek.r = z__1.r, ek.i = z__1.i;
L90:
    i__1 = kk;
    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)) 
	    != 0.) {
	i__2 = k;
	z_div(&z__1, &z__[k], &ap[kk]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
    }
    i__1 = kk;
    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)) 
	    == 0.) {
	i__2 = k;
	z__[i__2].r = 1., z__[i__2].i = 0.;
    }
    goto L110;
L100:
    km1k = ik + k - 1;
    km1km1 = ikm1 + k - 1;
    d_cnjg(&z__2, &ap[km1k]);
    z_div(&z__1, &ap[kk], &z__2);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &ap[km1km1], &ap[km1k]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    d_cnjg(&z__2, &ap[km1k]);
    z_div(&z__1, &z__[k], &z__2);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &z__[k - 1], &ap[km1k]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
L110:
    k -= ks;
    ik -= k;
    if (ks == 2) {
	ik -= k + 1;
    }
    goto L60;
L120:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*     SOLVE CTRANS(U)*Y = W */

    k = 1;
    ik = 0;
L130:
    if (k > *n) {
	goto L160;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L150;
    }
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotc_(&z__2, &i__3, &ap[ik + 1], &c__1, &z__[1], &c__1);
    z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    ikp1 = ik + k;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotc_(&z__2, &i__3, &ap[ikp1 + 1], &c__1, &z__[1], &c__1);
	z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
	z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L140;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L140:
L150:
    ik += k;
    if (ks == 2) {
	ik += k + 1;
    }
    k += ks;
    goto L130;
L160:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE U*D*V = Y */

    k = *n;
    ik = *n * (*n - 1) / 2;
L170:
    if (k == 0) {
	goto L230;
    }
    kk = ik + k;
    ikm1 = ik - (k - 1);
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == ks) {
	goto L190;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L180;
    }
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L180:
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &ap[ik + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - ks;
	caxpy_(&i__1, &z__[k - 1], &ap[ikm1 + 1], &c__1, &z__[1], &c__1);
    }
L190:
    if (ks == 2) {
	goto L210;
    }
    i__1 = k;
    i__2 = kk;
    if ((d__1 = z__[i__1].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(d__2)) 
	    <= (d__3 = ap[i__2].r, abs(d__3)) + (d__4 = d_imag(&ap[kk]), abs(
	    d__4))) {
	goto L200;
    }
    i__1 = kk;
    i__2 = k;
    s = ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)))
	     / ((d__3 = z__[i__2].r, abs(d__3)) + (d__4 = d_imag(&z__[k]), 
	    abs(d__4)));
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;
L200:
    i__1 = kk;
    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)) 
	    != 0.) {
	i__2 = k;
	z_div(&z__1, &z__[k], &ap[kk]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
    }
    i__1 = kk;
    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(d__2)) 
	    == 0.) {
	i__2 = k;
	z__[i__2].r = 1., z__[i__2].i = 0.;
    }
    goto L220;
L210:
    km1k = ik + k - 1;
    km1km1 = ikm1 + k - 1;
    d_cnjg(&z__2, &ap[km1k]);
    z_div(&z__1, &ap[kk], &z__2);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &ap[km1km1], &ap[km1k]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    d_cnjg(&z__2, &ap[km1k]);
    z_div(&z__1, &z__[k], &z__2);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &z__[k - 1], &ap[km1k]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
L220:
    k -= ks;
    ik -= k;
    if (ks == 2) {
	ik -= k + 1;
    }
    goto L170;
L230:
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     SOLVE CTRANS(U)*Z = V */

    k = 1;
    ik = 0;
L240:
    if (k > *n) {
	goto L270;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L260;
    }
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotc_(&z__2, &i__3, &ap[ik + 1], &c__1, &z__[1], &c__1);
    z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
    z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    ikp1 = ik + k;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotc_(&z__2, &i__3, &ap[ikp1 + 1], &c__1, &z__[1], &c__1);
	z__1.r = z__[i__2].r + z__2.r, z__1.i = z__[i__2].i + z__2.i;
	z__[i__1].r = z__1.r, z__[i__1].i = z__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L250;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L250:
L260:
    ik += k;
    if (ks == 2) {
	ik += k + 1;
    }
    k += ks;
    goto L240;
L270:
/*     MAKE ZNORM = 1.0 */
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* chpco_ */

/* Subroutine */ int chpfa_(doublecomplex *ap, integer *n, integer *kpvt, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal), d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k;
    static doublecomplex t, ak, bk;
    static integer ij, ik, jj, im, jk, kk, km1, km2, ijj, imj, imk;
    static doublecomplex akm1, bkm1;
    static integer ikm1, jkm1, km1k, imim, jmim, imax, jmax;
    static doublecomplex mulk;
    static logical swap;
    static doublereal alpha;
    static integer km1km1;
    static doublecomplex denom;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer kstep, imaxp1;
    static doublecomplex mulkm1;
    static doublereal absakk;
    extern integer icamax_(integer *, doublecomplex *, integer *);
    static doublereal colmax, rowmax;


/*     CHPFA FACTORS A COMPLEX HERMITIAN MATRIX STORED IN */
/*     PACKED FORM BY ELIMINATION WITH SYMMETRIC PIVOTING. */

/*     TO SOLVE  A*X = B , FOLLOW CHPFA BY CHPSL. */
/*     TO COMPUTE  INVERSE(A)*C , FOLLOW CHPFA BY CHPSL. */
/*     TO COMPUTE  DETERMINANT(A) , FOLLOW CHPFA BY CHPDI. */
/*     TO COMPUTE  INERTIA(A) , FOLLOW CHPFA BY CHPDI. */
/*     TO COMPUTE  INVERSE(A) , FOLLOW CHPFA BY CHPDI. */

/*     ON ENTRY */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                THE PACKED FORM OF A HERMITIAN MATRIX  A .  THE */
/*                COLUMNS OF THE UPPER TRIANGLE ARE STORED SEQUENTIALLY */
/*                IN A ONE-DIMENSIONAL ARRAY OF LENGTH  N*(N+1)/2 . */
/*                SEE COMMENTS BELOW FOR DETAILS. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*     OUTPUT */

/*        AP      A BLOCK DIAGONAL MATRIX AND THE MULTIPLIERS WHICH */
/*                WERE USED TO OBTAIN IT STORED IN PACKED FORM. */
/*                THE FACTORIZATION CAN BE WRITTEN  A = U*D*CTRANS(U) */
/*                WHERE  U  IS A PRODUCT OF PERMUTATION AND UNIT */
/*                UPPER TRIANGULAR MATRICES , CTRANS(U) IS THE */
/*                CONJUGATE TRANSPOSE OF  U , AND  D  IS BLOCK DIAGONAL */
/*                WITH 1 BY 1 AND 2 BY 2 BLOCKS. */

/*        KPVT    INTEGER(N) */
/*                AN INTEGER VECTOR OF PIVOT INDICES. */

/*        INFO    INTEGER */
/*                = 0  NORMAL VALUE. */
/*                = K  IF THE K-TH PIVOT BLOCK IS SINGULAR. THIS IS */
/*                     NOT AN ERROR CONDITION FOR THIS SUBROUTINE, */
/*                     BUT IT DOES INDICATE THAT CHPSL OR CHPDI MAY */
/*                     DIVIDE BY ZERO IF CALLED. */

/*     PACKED STORAGE */

/*          THE FOLLOWING PROGRAM SEGMENT WILL PACK THE UPPER */
/*          TRIANGLE OF A HERMITIAN MATRIX. */

/*                K = 0 */
/*                DO 20 J = 1, N */
/*                   DO 10 I = 1, J */
/*                      K = K + 1 */
/*                      AP(K)  = A(I,J) */
/*             10    CONTINUE */
/*             20 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CSWAP,ICAMAX */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,CONJG,DBLE,SQRT */

/*     INTERNAL VARIABLES */



/*     INITIALIZE */

/*     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE. */
    /* Parameter adjustments */
    --kpvt;
    --ap;

    /* Function Body */
    alpha = (sqrt(17.) + 1.) / 8.;

    *info = 0;

/*     MAIN LOOP ON K, WHICH GOES FROM N TO 1. */

    k = *n;
    ik = *n * (*n - 1) / 2;
L10:

/*        LEAVE THE LOOP IF K=0 OR K=1. */

/*     ...EXIT */
    if (k == 0) {
	goto L200;
    }
    if (k > 1) {
	goto L20;
    }
    kpvt[1] = 1;
    if ((d__1 = ap[1].r, abs(d__1)) + (d__2 = d_imag(&ap[1]), abs(d__2)) == 
	    0.) {
	*info = 1;
    }
/*     ......EXIT */
    goto L200;
L20:

/*        THIS SECTION OF CODE DETERMINES THE KIND OF */
/*        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED, */
/*        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND */
/*        SWAP WILL BE SET TO .TRUE. IF AN INTERCHANGE IS */
/*        REQUIRED. */

    km1 = k - 1;
    kk = ik + k;
    i__1 = kk;
    absakk = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[kk]), abs(
	    d__2));

/*        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*        COLUMN K. */

    i__1 = k - 1;
    imax = icamax_(&i__1, &ap[ik + 1], &c__1);
    imk = ik + imax;
    i__1 = imk;
    colmax = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[imk]), abs(
	    d__2));
    if (absakk < alpha * colmax) {
	goto L30;
    }
    kstep = 1;
    swap = FALSE_;
    goto L90;
L30:

/*           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN */
/*           ROW IMAX. */

    rowmax = 0.;
    imaxp1 = imax + 1;
    im = imax * (imax - 1) / 2;
    imj = im + (imax << 1);
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = imj;
	d__3 = rowmax, d__4 = (d__1 = ap[i__2].r, abs(d__1)) + (d__2 = d_imag(
		&ap[imj]), abs(d__2));
	rowmax = max(d__3,d__4);
	imj += j;
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = icamax_(&i__1, &ap[im + 1], &c__1);
    jmim = jmax + im;
/* Computing MAX */
    i__1 = jmim;
    d__3 = rowmax, d__4 = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[
	    jmim]), abs(d__2));
    rowmax = max(d__3,d__4);
L50:
    imim = imax + im;
    i__1 = imim;
    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = d_imag(&ap[imim]), abs(d__2))
	     < alpha * rowmax) {
	goto L60;
    }
    kstep = 1;
    swap = TRUE_;
    goto L80;
L60:
    if (absakk < alpha * colmax * (colmax / rowmax)) {
	goto L70;
    }
    kstep = 1;
    swap = FALSE_;
    goto L80;
L70:
    kstep = 2;
    swap = imax != km1;
L80:
L90:
    if (max(absakk,colmax) != 0.) {
	goto L100;
    }

/*           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP. */

    kpvt[k] = k;
    *info = k;
    goto L190;
L100:
    if (kstep == 2) {
	goto L140;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (! swap) {
	goto L120;
    }

/*              PERFORM AN INTERCHANGE. */

    cswap_(&imax, &ap[im + 1], &c__1, &ap[ik + 1], &c__1);
    imj = ik + imax;
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	jk = ik + j;
	d_cnjg(&z__1, &ap[jk]);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = jk;
	d_cnjg(&z__1, &ap[imj]);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = imj;
	ap[i__2].r = t.r, ap[i__2].i = t.i;
	imj -= j - 1;
/* L110: */
    }
L120:

/*           PERFORM THE ELIMINATION. */

    ij = ik - (k - 1);
    i__1 = km1;
    for (jj = 1; jj <= i__1; ++jj) {
	j = k - jj;
	jk = ik + j;
	i__2 = jk;
	z__2.r = -ap[i__2].r, z__2.i = -ap[i__2].i;
	z_div(&z__1, &z__2, &ap[kk]);
	mulk.r = z__1.r, mulk.i = z__1.i;
	d_cnjg(&z__1, &mulk);
	t.r = z__1.r, t.i = z__1.i;
	caxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	ijj = ij + j;
	i__2 = ijj;
	i__3 = ijj;
	d__1 = ap[i__3].r;
	z__1.r = d__1, z__1.i = 0.;
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = jk;
	ap[i__2].r = mulk.r, ap[i__2].i = mulk.i;
	ij -= j - 1;
/* L130: */
    }

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = k;
    if (swap) {
	kpvt[k] = imax;
    }
    goto L190;
L140:

/*           2 X 2 PIVOT BLOCK. */

    km1k = ik + k - 1;
    ikm1 = ik - (k - 1);
    if (! swap) {
	goto L160;
    }

/*              PERFORM AN INTERCHANGE. */

    cswap_(&imax, &ap[im + 1], &c__1, &ap[ikm1 + 1], &c__1);
    imj = ikm1 + imax;
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	jkm1 = ikm1 + j;
	d_cnjg(&z__1, &ap[jkm1]);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = jkm1;
	d_cnjg(&z__1, &ap[imj]);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = imj;
	ap[i__2].r = t.r, ap[i__2].i = t.i;
	imj -= j - 1;
/* L150: */
    }
    i__1 = km1k;
    t.r = ap[i__1].r, t.i = ap[i__1].i;
    i__1 = km1k;
    i__2 = imk;
    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
    i__1 = imk;
    ap[i__1].r = t.r, ap[i__1].i = t.i;
L160:

/*           PERFORM THE ELIMINATION. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    z_div(&z__1, &ap[kk], &ap[km1k]);
    ak.r = z__1.r, ak.i = z__1.i;
    km1km1 = ikm1 + k - 1;
    d_cnjg(&z__2, &ap[km1k]);
    z_div(&z__1, &ap[km1km1], &z__2);
    akm1.r = z__1.r, akm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = 1. - z__2.r, z__1.i = -z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    ij = ik - (k - 1) - (k - 2);
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	jk = ik + j;
	z_div(&z__1, &ap[jk], &ap[km1k]);
	bk.r = z__1.r, bk.i = z__1.i;
	jkm1 = ikm1 + j;
	d_cnjg(&z__2, &ap[km1k]);
	z_div(&z__1, &ap[jkm1], &z__2);
	bkm1.r = z__1.r, bkm1.i = z__1.i;
	z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + 
		akm1.i * bk.r;
	z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
	z_div(&z__1, &z__2, &denom);
	mulk.r = z__1.r, mulk.i = z__1.i;
	z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i 
		* bkm1.r;
	z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
	z_div(&z__1, &z__2, &denom);
	mulkm1.r = z__1.r, mulkm1.i = z__1.i;
	d_cnjg(&z__1, &mulk);
	t.r = z__1.r, t.i = z__1.i;
	caxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	d_cnjg(&z__1, &mulkm1);
	t.r = z__1.r, t.i = z__1.i;
	caxpy_(&j, &t, &ap[ikm1 + 1], &c__1, &ap[ij + 1], &c__1);
	i__2 = jk;
	ap[i__2].r = mulk.r, ap[i__2].i = mulk.i;
	i__2 = jkm1;
	ap[i__2].r = mulkm1.r, ap[i__2].i = mulkm1.i;
	ijj = ij + j;
	i__2 = ijj;
	i__3 = ijj;
	d__1 = ap[i__3].r;
	z__1.r = d__1, z__1.i = 0.;
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	ij -= j - 1;
/* L170: */
    }
L180:

/*           SET THE PIVOT ARRAY. */

    kpvt[k] = 1 - k;
    if (swap) {
	kpvt[k] = -imax;
    }
    kpvt[k - 1] = kpvt[k];
L190:
    ik -= k - 1;
    if (kstep == 2) {
	ik -= k - 2;
    }
    k -= kstep;
    goto L10;
L200:
    return 0;
} /* chpfa_ */

/* Subroutine */ int chpsl_(doublecomplex *ap, integer *n, integer *kpvt, 
	doublecomplex *b)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex ak, bk;
    static integer ik, kk, kp;
    static doublecomplex akm1, bkm1;
    static integer ikm1, km1k, ikp1;
    static doublecomplex temp;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static integer km1km1;
    static doublecomplex denom;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CHISL SOLVES THE COMPLEX HERMITIAN SYSTEM */
/*     A * X = B */
/*     USING THE FACTORS COMPUTED BY CHPFA. */

/*     ON ENTRY */

/*        AP      COMPLEX(N*(N+1)/2) */
/*                THE OUTPUT FROM CHPFA. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX  A . */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CHPFA. */

/*        B       COMPLEX(N) */
/*                THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       THE SOLUTION VECTOR  X . */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO MAY OCCUR IF  CHPCO  HAS SET RCOND .EQ. 0.0 */
/*        OR  CHPFA  HAS SET INFO .NE. 0  . */

/*     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX */
/*     WITH  P  COLUMNS */
/*           CALL CHPFA(AP,N,KPVT,INFO) */
/*           IF (INFO .NE. 0) GO TO ... */
/*           DO 10 J = 1, P */
/*              CALL CHPSL(AP,N,KPVT,C(1,J)) */
/*        10 CONTINUE */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CDOTC */
/*     FORTRAN CONJG,IABS */

/*     INTERNAL VARIABLES. */


/*     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND */
/*     D INVERSE TO B. */

    /* Parameter adjustments */
    --b;
    --kpvt;
    --ap;

    /* Function Body */
    k = *n;
    ik = *n * (*n - 1) / 2;
L10:
    if (k == 0) {
	goto L80;
    }
    kk = ik + k;
    if (kpvt[k] < 0) {
	goto L40;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L30;
    }
    kp = kpvt[k];
    if (kp == k) {
	goto L20;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L20:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 1;
    caxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
L30:

/*           APPLY D INVERSE. */

    i__1 = k;
    z_div(&z__1, &b[k], &ap[kk]);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    --k;
    ik -= k;
    goto L70;
L40:

/*           2 X 2 PIVOT BLOCK. */

    ikm1 = ik - (k - 1);
    if (k == 2) {
	goto L60;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k - 1) {
	goto L50;
    }

/*                 INTERCHANGE. */

    i__1 = k - 1;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k - 1;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L50:

/*              APPLY THE TRANSFORMATION. */

    i__1 = k - 2;
    caxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    caxpy_(&i__1, &b[k - 1], &ap[ikm1 + 1], &c__1, &b[1], &c__1);
L60:

/*           APPLY D INVERSE. */

    km1k = ik + k - 1;
    kk = ik + k;
    d_cnjg(&z__2, &ap[km1k]);
    z_div(&z__1, &ap[kk], &z__2);
    ak.r = z__1.r, ak.i = z__1.i;
    km1km1 = ikm1 + k - 1;
    z_div(&z__1, &ap[km1km1], &ap[km1k]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    d_cnjg(&z__2, &ap[km1k]);
    z_div(&z__1, &b[k], &z__2);
    bk.r = z__1.r, bk.i = z__1.i;
    z_div(&z__1, &b[k - 1], &ap[km1k]);
    bkm1.r = z__1.r, bkm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = z__2.r - 1., z__1.i = z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    i__1 = k;
    z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
    z_div(&z__1, &z__2, &denom);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    i__1 = k - 1;
    z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
    z_div(&z__1, &z__2, &denom);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    k += -2;
    ik = ik - (k + 1) - k;
L70:
    goto L10;
L80:

/*     LOOP FORWARD APPLYING THE TRANSFORMATIONS. */

    k = 1;
    ik = 0;
L90:
    if (k > *n) {
	goto L160;
    }
    if (kpvt[k] < 0) {
	goto L120;
    }

/*           1 X 1 PIVOT BLOCK. */

    if (k == 1) {
	goto L110;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotc_(&z__2, &i__3, &ap[ik + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    kp = kpvt[k];
    if (kp == k) {
	goto L100;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L100:
L110:
    ik += k;
    ++k;
    goto L150;
L120:

/*           2 X 2 PIVOT BLOCK. */

    if (k == 1) {
	goto L140;
    }

/*              APPLY THE TRANSFORMATION. */

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotc_(&z__2, &i__3, &ap[ik + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    ikp1 = ik + k;
    i__1 = k + 1;
    i__2 = k + 1;
    i__3 = k - 1;
    cdotc_(&z__2, &i__3, &ap[ikp1 + 1], &c__1, &b[1], &c__1);
    z__1.r = b[i__2].r + z__2.r, z__1.i = b[i__2].i + z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L130;
    }

/*                 INTERCHANGE. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L130:
L140:
    ik = ik + k + k + 1;
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* chpsl_ */

/* Subroutine */ int chpdi_(doublecomplex *ap, integer *n, integer *kpvt, 
	doublereal *det, integer *inert, doublecomplex *work, integer *job)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static doublereal d__;
    static integer j, k;
    static doublereal t, ak;
    static integer jb, ij, ik, jk, kk, ks, km1;
    static doublereal ten;
    static integer iks, ksj;
    static doublereal akp1;
    static integer ikp1, jkp1, kkp1;
    static doublecomplex temp, akkp1;
    static integer kskp1;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical nodet;
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), cswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), caxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);
    static integer kstep;
    static logical noert, noinv;


/*     CHPDI COMPUTES THE DETERMINANT, INERTIA AND INVERSE */
/*     OF A COMPLEX HERMITIAN MATRIX USING THE FACTORS FROM CHPFA, */
/*     WHERE THE MATRIX IS STORED IN PACKED FORM. */

/*     ON ENTRY */

/*        AP      COMPLEX (N*(N+1)/2) */
/*                THE OUTPUT FROM CHPFA. */

/*        N       INTEGER */
/*                THE ORDER OF THE MATRIX A. */

/*        KPVT    INTEGER(N) */
/*                THE PIVOT VECTOR FROM CHPFA. */

/*        WORK    COMPLEX(N) */
/*                WORK VECTOR.  CONTENTS IGNORED. */

/*        JOB     INTEGER */
/*                JOB HAS THE DECIMAL EXPANSION  ABC  WHERE */
/*                   IF  C .NE. 0, THE INVERSE IS COMPUTED, */
/*                   IF  B .NE. 0, THE DETERMINANT IS COMPUTED, */
/*                   IF  A .NE. 0, THE INERTIA IS COMPUTED. */

/*                FOR EXAMPLE, JOB = 111  GIVES ALL THREE. */

/*     ON RETURN */

/*        VARIABLES NOT REQUESTED BY JOB ARE NOT USED. */

/*        AP     CONTAINS THE UPPER TRIANGLE OF THE INVERSE OF */
/*               THE ORIGINAL MATRIX, STORED IN PACKED FORM. */
/*               THE COLUMNS OF THE UPPER TRIANGLE ARE STORED */
/*               SEQUENTIALLY IN A ONE-DIMENSIONAL ARRAY. */

/*        DET    REAL(2) */
/*               DETERMINANT OF ORIGINAL MATRIX. */
/*               DETERMINANT = DET(1) * 10.0**DET(2) */
/*               WITH 1.0 .LE. ABS(DET(1)) .LT. 10.0 */
/*               OR DET(1) = 0.0. */

/*        INERT  INTEGER(3) */
/*               THE INERTIA OF THE ORIGINAL MATRIX. */
/*               INERT(1)  =  NUMBER OF POSITIVE EIGENVALUES. */
/*               INERT(2)  =  NUMBER OF NEGATIVE EIGENVALUES. */
/*               INERT(3)  =  NUMBER OF ZERO EIGENVALUES. */

/*     ERROR CONDITION */

/*        A DIVISION BY ZERO WILL OCCUR IF THE INVERSE IS REQUESTED */
/*        AND  CHPCO  HAS SET RCOND .EQ. 0.0 */
/*        OR  CHPFA  HAS SET  INFO .NE. 0 . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JAMES BUNCH, UNIV. CALIF. SAN DIEGO, ARGONNE NAT. LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CCOPY,CDOTC,CSWAP */
/*     FORTRAN ABS,CABS,DCMPLX,CONJG,IABS,MOD,DBLE */

/*     INTERNAL VARIABLES. */


    /* Parameter adjustments */
    --work;
    --inert;
    --det;
    --kpvt;
    --ap;

    /* Function Body */
    noinv = *job % 10 == 0;
    nodet = *job % 100 / 10 == 0;
    noert = *job % 1000 / 100 == 0;

    if (nodet && noert) {
	goto L140;
    }
    if (noert) {
	goto L10;
    }
    inert[1] = 0;
    inert[2] = 0;
    inert[3] = 0;
L10:
    if (nodet) {
	goto L20;
    }
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
L20:
    t = 0.;
    ik = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk = ik + k;
	i__2 = kk;
	d__ = ap[i__2].r;

/*           CHECK IF 1 BY 1 */

	if (kpvt[k] > 0) {
	    goto L50;
	}

/*              2 BY 2 BLOCK */
/*              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = ABS(S) */
/*                      (S  C) */
/*              TO AVOID UNDERFLOW/OVERFLOW TROUBLES. */
/*              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG. */

	if (t != 0.) {
	    goto L30;
	}
	ikp1 = ik + k;
	kkp1 = ikp1 + k;
	t = z_abs(&ap[kkp1]);
	i__2 = kkp1 + 1;
	d__ = d__ / t * ap[i__2].r - t;
	goto L40;
L30:
	d__ = t;
	t = 0.;
L40:
L50:

	if (noert) {
	    goto L60;
	}
	if (d__ > 0.) {
	    ++inert[1];
	}
	if (d__ < 0.) {
	    ++inert[2];
	}
	if (d__ == 0.) {
	    ++inert[3];
	}
L60:

	if (nodet) {
	    goto L120;
	}
	det[1] = d__ * det[1];
	if (det[1] == 0.) {
	    goto L110;
	}
L70:
	if (abs(det[1]) >= 1.) {
	    goto L80;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L70;
L80:
L90:
	if (abs(det[1]) < ten) {
	    goto L100;
	}
	det[1] /= ten;
	det[2] += 1.;
	goto L90;
L100:
L110:
L120:
	ik += k;
/* L130: */
    }
L140:

/*     COMPUTE INVERSE(A) */

    if (noinv) {
	goto L270;
    }
    k = 1;
    ik = 0;
L150:
    if (k > *n) {
	goto L260;
    }
    km1 = k - 1;
    kk = ik + k;
    ikp1 = ik + k;
    kkp1 = ikp1 + k;
    if (kpvt[k] < 0) {
	goto L180;
    }

/*              1 BY 1 */

    i__1 = kk;
    i__2 = kk;
    d__1 = 1. / ap[i__2].r;
    z__1.r = d__1, z__1.i = 0.;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L170;
    }
    ccopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	i__2 = jk;
	cdotc_(&z__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L160: */
    }
    i__1 = kk;
    i__2 = kk;
    cdotc_(&z__3, &km1, &work[1], &c__1, &ap[ik + 1], &c__1);
    d__1 = z__3.r;
    z__2.r = d__1, z__2.i = 0.;
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
L170:
    kstep = 1;
    goto L220;
L180:

/*              2 BY 2 */

    t = z_abs(&ap[kkp1]);
    i__1 = kk;
    ak = ap[i__1].r / t;
    i__1 = kkp1 + 1;
    akp1 = ap[i__1].r / t;
    i__1 = kkp1;
    z__1.r = ap[i__1].r / t, z__1.i = ap[i__1].i / t;
    akkp1.r = z__1.r, akkp1.i = z__1.i;
    d__ = t * (ak * akp1 - 1.);
    i__1 = kk;
    d__1 = akp1 / d__;
    z__1.r = d__1, z__1.i = 0.;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    i__1 = kkp1 + 1;
    d__1 = ak / d__;
    z__1.r = d__1, z__1.i = 0.;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    i__1 = kkp1;
    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
    z__1.r = z__2.r / d__, z__1.i = z__2.i / d__;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L210;
    }
    ccopy_(&km1, &ap[ikp1 + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jkp1 = ikp1 + j;
	i__2 = jkp1;
	cdotc_(&z__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ikp1 + 1], &c__1);
	ij += j;
/* L190: */
    }
    i__1 = kkp1 + 1;
    i__2 = kkp1 + 1;
    cdotc_(&z__3, &km1, &work[1], &c__1, &ap[ikp1 + 1], &c__1);
    d__1 = z__3.r;
    z__2.r = d__1, z__2.i = 0.;
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    i__1 = kkp1;
    i__2 = kkp1;
    cdotc_(&z__2, &km1, &ap[ik + 1], &c__1, &ap[ikp1 + 1], &c__1);
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    ccopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	i__2 = jk;
	cdotc_(&z__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L200: */
    }
    i__1 = kk;
    i__2 = kk;
    cdotc_(&z__3, &km1, &work[1], &c__1, &ap[ik + 1], &c__1);
    d__1 = z__3.r;
    z__2.r = d__1, z__2.i = 0.;
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
L210:
    kstep = 2;
L220:

/*           SWAP */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L250;
    }
    iks = ks * (ks - 1) / 2;
    cswap_(&ks, &ap[iks + 1], &c__1, &ap[ik + 1], &c__1);
    ksj = ik + ks;
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	jk = ik + j;
	d_cnjg(&z__1, &ap[jk]);
	temp.r = z__1.r, temp.i = z__1.i;
	i__2 = jk;
	d_cnjg(&z__1, &ap[ksj]);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = ksj;
	ap[i__2].r = temp.r, ap[i__2].i = temp.i;
	ksj -= j - 1;
/* L230: */
    }
    if (kstep == 1) {
	goto L240;
    }
    kskp1 = ikp1 + ks;
    i__1 = kskp1;
    temp.r = ap[i__1].r, temp.i = ap[i__1].i;
    i__1 = kskp1;
    i__2 = kkp1;
    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
    i__1 = kkp1;
    ap[i__1].r = temp.r, ap[i__1].i = temp.i;
L240:
L250:
    ik += k;
    if (kstep == 2) {
	ik = ik + k + 1;
    }
    k += kstep;
    goto L150;
L260:
L270:
    return 0;
} /* chpdi_ */

/* Subroutine */ int ctrco_(doublecomplex *t, integer *ldt, integer *n, 
	doublereal *rcond, doublecomplex *z__, integer *job)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, l;
    static doublereal s;
    static doublecomplex w;
    static integer i1, j1, j2;
    static doublecomplex ek;
    static integer kk;
    static doublereal sm;
    static doublecomplex wk, wkm;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical lower;
    static doublereal tnorm, ynorm;
    extern /* Subroutine */ int csscal_(integer *, doublereal *, 
	    doublecomplex *, integer *);
    extern doublereal scasum_(integer *, doublecomplex *, integer *);


/*     CTRCO ESTIMATES THE CONDITION OF A COMPLEX TRIANGULAR MATRIX. */

/*     ON ENTRY */

/*        T       COMPLEX(LDT,N) */
/*                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO */
/*                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND */
/*                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE */
/*                USED TO STORE OTHER INFORMATION. */

/*        LDT     INTEGER */
/*                LDT IS THE LEADING DIMENSION OF THE ARRAY T. */

/*        N       INTEGER */
/*                N IS THE ORDER OF THE SYSTEM. */

/*        JOB     INTEGER */
/*                = 0         T  IS LOWER TRIANGULAR. */
/*                = NONZERO   T  IS UPPER TRIANGULAR. */

/*     ON RETURN */

/*        RCOND   REAL */
/*                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  T . */
/*                FOR THE SYSTEM  T*X = B , RELATIVE PERTURBATIONS */
/*                IN  T  AND  B  OF SIZE  EPSILON  MAY CAUSE */
/*                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND . */
/*                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION */
/*                           1.0 + RCOND .EQ. 1.0 */
/*                IS TRUE, THEN  T  MAY BE SINGULAR TO WORKING */
/*                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF */
/*                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE */
/*                UNDERFLOWS. */

/*        Z       COMPLEX(N) */
/*                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT. */
/*                IF  T  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS */
/*                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT */
/*                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) . */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CSSCAL,SCASUM */
/*     FORTRAN ABS,DIMAG,MAX,DCMPLX,CONJG,DBLE */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --z__;

    /* Function Body */
    lower = *job == 0;

/*     COMPUTE 1-NORM OF T */

    tnorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = j;
	if (lower) {
	    l = *n + 1 - j;
	}
	i1 = 1;
	if (lower) {
	    i1 = j;
	}
/* Computing MAX */
	d__1 = tnorm, d__2 = scasum_(&l, &t[i1 + j * t_dim1], &c__1);
	tnorm = max(d__1,d__2);
/* L10: */
    }

/*     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) . */
/*     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  CTRANS(T)*Y = E . */
/*     CTRANS(T)  IS THE CONJUGATE TRANSPOSE OF T . */
/*     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL */
/*     GROWTH IN THE ELEMENTS OF Y . */
/*     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW. */

/*     SOLVE CTRANS(T)*Y = E */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L20: */
    }
    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = kk;
	if (lower) {
	    k = *n + 1 - kk;
	}
	i__2 = k;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) != 0.) {
	    i__3 = k;
	    z__2.r = -z__[i__3].r, z__2.i = -z__[i__3].i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = d_imag(&ek), abs(d__4));
	    d__8 = (d__5 = z__1.r, abs(d__5)) + (d__6 = d_imag(&z__1), abs(
		    d__6));
	    z__4.r = z__1.r / d__8, z__4.i = z__1.i / d__8;
	    z__3.r = d__7 * z__4.r, z__3.i = d__7 * z__4.i;
	    ek.r = z__3.r, ek.i = z__3.i;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = k + k * t_dim1;
	if ((d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(d__2)) <= 
		(d__3 = t[i__3].r, abs(d__3)) + (d__4 = d_imag(&t[k + k * 
		t_dim1]), abs(d__4))) {
	    goto L30;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = k + k * t_dim1;
	s = ((d__1 = t[i__3].r, abs(d__1)) + (d__2 = d_imag(&t[k + k * t_dim1]
		), abs(d__2))) / ((d__3 = z__1.r, abs(d__3)) + (d__4 = d_imag(
		&z__1), abs(d__4)));
	csscal_(n, &s, &z__[1], &c__1);
	z__2.r = s, z__2.i = 0.;
	z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + 
		z__2.i * ek.r;
	ek.r = z__1.r, ek.i = z__1.i;
L30:
	i__2 = k;
	z__1.r = ek.r - z__[i__2].r, z__1.i = ek.i - z__[i__2].i;
	wk.r = z__1.r, wk.i = z__1.i;
	z__2.r = -ek.r, z__2.i = -ek.i;
	i__2 = k;
	z__1.r = z__2.r - z__[i__2].r, z__1.i = z__2.i - z__[i__2].i;
	wkm.r = z__1.r, wkm.i = z__1.i;
	s = (d__1 = wk.r, abs(d__1)) + (d__2 = d_imag(&wk), abs(d__2));
	sm = (d__1 = wkm.r, abs(d__1)) + (d__2 = d_imag(&wkm), abs(d__2));
	i__2 = k + k * t_dim1;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = d_imag(&t[k + k * t_dim1])
		, abs(d__2)) == 0.) {
	    goto L40;
	}
	d_cnjg(&z__2, &t[k + k * t_dim1]);
	z_div(&z__1, &wk, &z__2);
	wk.r = z__1.r, wk.i = z__1.i;
	d_cnjg(&z__2, &t[k + k * t_dim1]);
	z_div(&z__1, &wkm, &z__2);
	wkm.r = z__1.r, wkm.i = z__1.i;
	goto L50;
L40:
	wk.r = 1., wk.i = 0.;
	wkm.r = 1., wkm.i = 0.;
L50:
	if (kk == *n) {
	    goto L90;
	}
	j1 = k + 1;
	if (lower) {
	    j1 = 1;
	}
	j2 = *n;
	if (lower) {
	    j2 = k - 1;
	}
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    i__3 = j;
	    d_cnjg(&z__4, &t[k + j * t_dim1]);
	    z__3.r = wkm.r * z__4.r - wkm.i * z__4.i, z__3.i = wkm.r * z__4.i 
		    + wkm.i * z__4.r;
	    z__2.r = z__[i__3].r + z__3.r, z__2.i = z__[i__3].i + z__3.i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    sm += (d__1 = z__1.r, abs(d__1)) + (d__2 = d_imag(&z__1), abs(
		    d__2));
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &t[k + j * t_dim1]);
	    z__2.r = wk.r * z__3.r - wk.i * z__3.i, z__2.i = wk.r * z__3.i + 
		    wk.i * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = j;
	    s += (d__1 = z__[i__3].r, abs(d__1)) + (d__2 = d_imag(&z__[j]), 
		    abs(d__2));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	z__1.r = wkm.r - wk.r, z__1.i = wkm.i - wk.i;
	w.r = z__1.r, w.i = z__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &t[k + j * t_dim1]);
	    z__2.r = w.r * z__3.r - w.i * z__3.i, z__2.i = w.r * z__3.i + w.i 
		    * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
/* L70: */
	}
L80:
L90:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L100: */
    }
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     SOLVE T*Z = Y */

    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = *n + 1 - kk;
	if (lower) {
	    k = kk;
	}
	i__2 = k;
	i__3 = k + k * t_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = d_imag(&z__[k]), abs(
		d__2)) <= (d__3 = t[i__3].r, abs(d__3)) + (d__4 = d_imag(&t[k 
		+ k * t_dim1]), abs(d__4))) {
	    goto L110;
	}
	i__2 = k + k * t_dim1;
	i__3 = k;
	s = ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = d_imag(&t[k + k * t_dim1]
		), abs(d__2))) / ((d__3 = z__[i__3].r, abs(d__3)) + (d__4 = 
		d_imag(&z__[k]), abs(d__4)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L110:
	i__2 = k + k * t_dim1;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = d_imag(&t[k + k * t_dim1])
		, abs(d__2)) != 0.) {
	    i__3 = k;
	    z_div(&z__1, &z__[k], &t[k + k * t_dim1]);
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	}
	i__2 = k + k * t_dim1;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = d_imag(&t[k + k * t_dim1])
		, abs(d__2)) == 0.) {
	    i__3 = k;
	    z__[i__3].r = 1., z__[i__3].i = 0.;
	}
	i1 = 1;
	if (lower) {
	    i1 = k + 1;
	}
	if (kk >= *n) {
	    goto L120;
	}
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	w.r = z__1.r, w.i = z__1.i;
	i__2 = *n - kk;
	caxpy_(&i__2, &w, &t[i1 + k * t_dim1], &c__1, &z__[i1], &c__1);
L120:
/* L130: */
	;
    }
/*     MAKE ZNORM = 1.0 */
    s = 1. / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (tnorm != 0.) {
	*rcond = ynorm / tnorm;
    }
    if (tnorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* ctrco_ */

/* Subroutine */ int ctrsl_(doublecomplex *t, integer *ldt, integer *n, 
	doublecomplex *b, integer *job, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, jj, case__;
    static doublecomplex temp;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);



/*     CTRSL SOLVES SYSTEMS OF THE FORM */

/*                   T * X = B */
/*     OR */
/*                   CTRANS(T) * X = B */

/*     WHERE T IS A TRIANGULAR MATRIX OF ORDER N. HERE CTRANS(T) */
/*     DENOTES THE CONJUGATE TRANSPOSE OF THE MATRIX T. */

/*     ON ENTRY */

/*         T         COMPLEX(LDT,N) */
/*                   T CONTAINS THE MATRIX OF THE SYSTEM. THE ZERO */
/*                   ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND */
/*                   THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE */
/*                   USED TO STORE OTHER INFORMATION. */

/*         LDT       INTEGER */
/*                   LDT IS THE LEADING DIMENSION OF THE ARRAY T. */

/*         N         INTEGER */
/*                   N IS THE ORDER OF THE SYSTEM. */

/*         B         COMPLEX(N). */
/*                   B CONTAINS THE RIGHT HAND SIDE OF THE SYSTEM. */

/*         JOB       INTEGER */
/*                   JOB SPECIFIES WHAT KIND OF SYSTEM IS TO BE SOLVED. */
/*                   IF JOB IS */

/*                        00   SOLVE T*X=B, T LOWER TRIANGULAR, */
/*                        01   SOLVE T*X=B, T UPPER TRIANGULAR, */
/*                        10   SOLVE CTRANS(T)*X=B, T LOWER TRIANGULAR, */
/*                        11   SOLVE CTRANS(T)*X=B, T UPPER TRIANGULAR. */

/*     ON RETURN */

/*         B         B CONTAINS THE SOLUTION, IF INFO .EQ. 0. */
/*                   OTHERWISE B IS UNALTERED. */

/*         INFO      INTEGER */
/*                   INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR. */
/*                   OTHERWISE INFO CONTAINS THE INDEX OF */
/*                   THE FIRST ZERO DIAGONAL ELEMENT OF T. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G. W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CDOTC */
/*     FORTRAN ABS,DIMAG,CONJG,MOD,DBLE */

/*     INTERNAL VARIABLES */


/*     BEGIN BLOCK PERMITTING ...EXITS TO 150 */

/*        CHECK FOR ZERO DIAGONAL ELEMENTS. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (*info = 1; *info <= i__1; ++(*info)) {
/*     ......EXIT */
	i__2 = *info + *info * t_dim1;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = d_imag(&t[*info + *info * 
		t_dim1]), abs(d__2)) == 0.) {
	    goto L150;
	}
/* L10: */
    }
    *info = 0;

/*        DETERMINE THE TASK AND GO TO IT. */

    case__ = 1;
    if (*job % 10 != 0) {
	case__ = 2;
    }
    if (*job % 100 / 10 != 0) {
	case__ += 2;
    }
    switch (case__) {
	case 1:  goto L20;
	case 2:  goto L50;
	case 3:  goto L80;
	case 4:  goto L110;
    }

/*        SOLVE T*X=B FOR T LOWER TRIANGULAR */

L20:
    z_div(&z__1, &b[1], &t[t_dim1 + 1]);
    b[1].r = z__1.r, b[1].i = z__1.i;
    if (*n < 2) {
	goto L40;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	temp.r = z__1.r, temp.i = z__1.i;
	i__2 = *n - j + 1;
	caxpy_(&i__2, &temp, &t[j + (j - 1) * t_dim1], &c__1, &b[j], &c__1);
	i__2 = j;
	z_div(&z__1, &b[j], &t[j + j * t_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L30: */
    }
L40:
    goto L140;

/*        SOLVE T*X=B FOR T UPPER TRIANGULAR. */

L50:
    i__1 = *n;
    z_div(&z__1, &b[*n], &t[*n + *n * t_dim1]);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    if (*n < 2) {
	goto L70;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	i__2 = j + 1;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	temp.r = z__1.r, temp.i = z__1.i;
	caxpy_(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c__1, &b[1], &c__1);
	i__2 = j;
	z_div(&z__1, &b[j], &t[j + j * t_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L60: */
    }
L70:
    goto L140;

/*        SOLVE CTRANS(T)*X=B FOR T LOWER TRIANGULAR. */

L80:
    i__1 = *n;
    d_cnjg(&z__2, &t[*n + *n * t_dim1]);
    z_div(&z__1, &b[*n], &z__2);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    if (*n < 2) {
	goto L100;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	i__2 = j;
	i__3 = j;
	i__4 = jj - 1;
	cdotc_(&z__2, &i__4, &t[j + 1 + j * t_dim1], &c__1, &b[j + 1], &c__1);
	z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	i__2 = j;
	d_cnjg(&z__2, &t[j + j * t_dim1]);
	z_div(&z__1, &b[j], &z__2);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L90: */
    }
L100:
    goto L140;

/*        SOLVE CTRANS(T)*X=B FOR T UPPER TRIANGULAR. */

L110:
    d_cnjg(&z__2, &t[t_dim1 + 1]);
    z_div(&z__1, &b[1], &z__2);
    b[1].r = z__1.r, b[1].i = z__1.i;
    if (*n < 2) {
	goto L130;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	i__4 = j - 1;
	cdotc_(&z__2, &i__4, &t[j * t_dim1 + 1], &c__1, &b[1], &c__1);
	z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	i__2 = j;
	d_cnjg(&z__2, &t[j + j * t_dim1]);
	z_div(&z__1, &b[j], &z__2);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L120: */
    }
L130:
L140:
L150:
    return 0;
} /* ctrsl_ */

/* Subroutine */ int ctrdi_(doublecomplex *t, integer *ldt, integer *n, 
	doublecomplex *det, integer *job, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j, k, kb, km1, kp1;
    static doublereal ten;
    static doublecomplex temp;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CTRDI COMPUTES THE DETERMINANT AND INVERSE OF A COMPLEX */
/*     TRIANGULAR MATRIX. */

/*     ON ENTRY */

/*        T       COMPLEX(LDT,N) */
/*                T CONTAINS THE TRIANGULAR MATRIX. THE ZERO */
/*                ELEMENTS OF THE MATRIX ARE NOT REFERENCED, AND */
/*                THE CORRESPONDING ELEMENTS OF THE ARRAY CAN BE */
/*                USED TO STORE OTHER INFORMATION. */

/*        LDT     INTEGER */
/*                LDT IS THE LEADING DIMENSION OF THE ARRAY T. */

/*        N       INTEGER */
/*                N IS THE ORDER OF THE SYSTEM. */

/*        JOB     INTEGER */
/*                = 010       NO DET, INVERSE OF LOWER TRIANGULAR. */
/*                = 011       NO DET, INVERSE OF UPPER TRIANGULAR. */
/*                = 100       DET, NO INVERSE. */
/*                = 110       DET, INVERSE OF LOWER TRIANGULAR. */
/*                = 111       DET, INVERSE OF UPPER TRIANGULAR. */

/*     ON RETURN */

/*        T       INVERSE OF ORIGINAL MATRIX IF REQUESTED. */
/*                OTHERWISE UNCHANGED. */

/*        DET     COMPLEX(2) */
/*                DETERMINANT OF ORIGINAL MATRIX IF REQUESTED. */
/*                OTHERWISE NOT REFERENCED. */
/*                DETERMINANT = DET(1) * 10.0**DET(2) */
/*                WITH  1.0 .LE. CABS1(DET(1)) .LT. 10.0 */
/*                OR  DET(1) .EQ. 0.0 . */

/*        INFO    INTEGER */
/*                INFO CONTAINS ZERO IF THE SYSTEM IS NONSINGULAR */
/*                AND THE INVERSE IS REQUESTED. */
/*                OTHERWISE INFO CONTAINS THE INDEX OF */
/*                A ZERO DIAGONAL ELEMENT OF T. */


/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB. */

/*     SUBROUTINES AND FUNCTIONS */

/*     BLAS CAXPY,CSCAL */
/*     FORTRAN ABS,DIMAG,DCMPLX,MOD,DBLE */

/*     INTERNAL VARIABLES */


/*     BEGIN BLOCK PERMITTING ...EXITS TO 180 */

/*        COMPUTE DETERMINANT */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --det;

    /* Function Body */
    if (*job / 100 == 0) {
	goto L70;
    }
    det[1].r = 1., det[1].i = 0.;
    det[2].r = 0., det[2].i = 0.;
    ten = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + i__ * t_dim1;
	z__1.r = t[i__2].r * det[1].r - t[i__2].i * det[1].i, z__1.i = t[i__2]
		.r * det[1].i + t[i__2].i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
/*           ...EXIT */
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 == 0.) {
	    goto L60;
	}
L10:
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 >= 1.) {
	    goto L20;
	}
	z__2.r = ten, z__2.i = 0.;
	z__1.r = z__2.r * det[1].r - z__2.i * det[1].i, z__1.i = z__2.r * det[
		1].i + z__2.i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r - 1., z__1.i = det[2].i - 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L10;
L20:
L30:
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = d_imag(&det[1]), abs(d__2))
		 < ten) {
	    goto L40;
	}
	z__2.r = ten, z__2.i = 0.;
	z_div(&z__1, &det[1], &z__2);
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r + 1., z__1.i = det[2].i + 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*        COMPUTE INVERSE OF UPPER TRIANGULAR */

    if (*job / 10 % 10 == 0) {
	goto L170;
    }
    if (*job % 10 == 0) {
	goto L120;
    }
/*              BEGIN BLOCK PERMITTING ...EXITS TO 110 */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	*info = k;
/*              ......EXIT */
	i__2 = k + k * t_dim1;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = d_imag(&t[k + k * t_dim1])
		, abs(d__2)) == 0.) {
	    goto L110;
	}
	i__2 = k + k * t_dim1;
	z_div(&z__1, &c_b1092, &t[k + k * t_dim1]);
	t[i__2].r = z__1.r, t[i__2].i = z__1.i;
	i__2 = k + k * t_dim1;
	z__1.r = -t[i__2].r, z__1.i = -t[i__2].i;
	temp.r = z__1.r, temp.i = z__1.i;
	i__2 = k - 1;
	cscal_(&i__2, &temp, &t[k * t_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = k + j * t_dim1;
	    temp.r = t[i__3].r, temp.i = t[i__3].i;
	    i__3 = k + j * t_dim1;
	    t[i__3].r = 0., t[i__3].i = 0.;
	    caxpy_(&k, &temp, &t[k * t_dim1 + 1], &c__1, &t[j * t_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }
    *info = 0;
L110:
    goto L160;
L120:

/*              COMPUTE INVERSE OF LOWER TRIANGULAR */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	*info = k;
/*     ............EXIT */
	i__2 = k + k * t_dim1;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = d_imag(&t[k + k * t_dim1])
		, abs(d__2)) == 0.) {
	    goto L180;
	}
	i__2 = k + k * t_dim1;
	z_div(&z__1, &c_b1092, &t[k + k * t_dim1]);
	t[i__2].r = z__1.r, t[i__2].i = z__1.i;
	i__2 = k + k * t_dim1;
	z__1.r = -t[i__2].r, z__1.i = -t[i__2].i;
	temp.r = z__1.r, temp.i = z__1.i;
	if (k != *n) {
	    i__2 = *n - k;
	    cscal_(&i__2, &temp, &t[k + 1 + k * t_dim1], &c__1);
	}
	km1 = k - 1;
	if (km1 < 1) {
	    goto L140;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = k + j * t_dim1;
	    temp.r = t[i__3].r, temp.i = t[i__3].i;
	    i__3 = k + j * t_dim1;
	    t[i__3].r = 0., t[i__3].i = 0.;
	    i__3 = *n - k + 1;
	    caxpy_(&i__3, &temp, &t[k + k * t_dim1], &c__1, &t[k + j * t_dim1]
		    , &c__1);
/* L130: */
	}
L140:
/* L150: */
	;
    }
    *info = 0;
L160:
L170:
L180:
    return 0;
} /* ctrdi_ */

/* Subroutine */ int cgtsl_(integer *n, doublecomplex *c__, doublecomplex *
	d__, doublecomplex *e, doublecomplex *b, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3, z__4, z__5;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex t;
    static integer kb, kp1, nm1, nm2;


/*     CGTSL GIVEN A GENERAL TRIDIAGONAL MATRIX AND A RIGHT HAND */
/*     SIDE WILL FIND THE SOLUTION. */

/*     ON ENTRY */

/*        N       INTEGER */
/*                IS THE ORDER OF THE TRIDIAGONAL MATRIX. */

/*        C       COMPLEX(N) */
/*                IS THE SUBDIAGONAL OF THE TRIDIAGONAL MATRIX. */
/*                C(2) THROUGH C(N) SHOULD CONTAIN THE SUBDIAGONAL. */
/*                ON OUTPUT C IS DESTROYED. */

/*        D       COMPLEX(N) */
/*                IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX. */
/*                ON OUTPUT D IS DESTROYED. */

/*        E       COMPLEX(N) */
/*                IS THE SUPERDIAGONAL OF THE TRIDIAGONAL MATRIX. */
/*                E(1) THROUGH E(N-1) SHOULD CONTAIN THE SUPERDIAGONAL. */
/*                ON OUTPUT E IS DESTROYED. */

/*        B       COMPLEX(N) */
/*                IS THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B       IS THE SOLUTION VECTOR. */

/*        INFO    INTEGER */
/*                = 0 NORMAL VALUE. */
/*                = K IF THE K-TH ELEMENT OF THE DIAGONAL BECOMES */
/*                    EXACTLY ZERO.  THE SUBROUTINE RETURNS WHEN */
/*                    THIS IS DETECTED. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JACK DONGARRA, ARGONNE NATIONAL LABORATORY. */

/*     NO EXTERNALS */
/*     FORTRAN ABS,DIMAG,DBLE */

/*     INTERNAL VARIABLES */

/*     BEGIN BLOCK PERMITTING ...EXITS TO 100 */

    /* Parameter adjustments */
    --b;
    --e;
    --d__;
    --c__;

    /* Function Body */
    *info = 0;
    c__[1].r = d__[1].r, c__[1].i = d__[1].i;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L40;
    }
    d__[1].r = e[1].r, d__[1].i = e[1].i;
    e[1].r = 0., e[1].i = 0.;
    i__1 = *n;
    e[i__1].r = 0., e[i__1].i = 0.;

    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*              FIND THE LARGEST OF THE TWO ROWS */

	i__2 = kp1;
	i__3 = k;
	if ((d__1 = c__[i__2].r, abs(d__1)) + (d__2 = d_imag(&c__[kp1]), abs(
		d__2)) < (d__3 = c__[i__3].r, abs(d__3)) + (d__4 = d_imag(&
		c__[k]), abs(d__4))) {
	    goto L10;
	}

/*                 INTERCHANGE ROW */

	i__2 = kp1;
	t.r = c__[i__2].r, t.i = c__[i__2].i;
	i__2 = kp1;
	i__3 = k;
	c__[i__2].r = c__[i__3].r, c__[i__2].i = c__[i__3].i;
	i__2 = k;
	c__[i__2].r = t.r, c__[i__2].i = t.i;
	i__2 = kp1;
	t.r = d__[i__2].r, t.i = d__[i__2].i;
	i__2 = kp1;
	i__3 = k;
	d__[i__2].r = d__[i__3].r, d__[i__2].i = d__[i__3].i;
	i__2 = k;
	d__[i__2].r = t.r, d__[i__2].i = t.i;
	i__2 = kp1;
	t.r = e[i__2].r, t.i = e[i__2].i;
	i__2 = kp1;
	i__3 = k;
	e[i__2].r = e[i__3].r, e[i__2].i = e[i__3].i;
	i__2 = k;
	e[i__2].r = t.r, e[i__2].i = t.i;
	i__2 = kp1;
	t.r = b[i__2].r, t.i = b[i__2].i;
	i__2 = kp1;
	i__3 = k;
	b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
	i__2 = k;
	b[i__2].r = t.r, b[i__2].i = t.i;
L10:

/*              ZERO ELEMENTS */

	i__2 = k;
	if ((d__1 = c__[i__2].r, abs(d__1)) + (d__2 = d_imag(&c__[k]), abs(
		d__2)) != 0.) {
	    goto L20;
	}
	*info = k;
/*     ............EXIT */
	goto L100;
L20:
	i__2 = kp1;
	z__2.r = -c__[i__2].r, z__2.i = -c__[i__2].i;
	z_div(&z__1, &z__2, &c__[k]);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = kp1;
	i__3 = kp1;
	i__4 = k;
	z__2.r = t.r * d__[i__4].r - t.i * d__[i__4].i, z__2.i = t.r * d__[
		i__4].i + t.i * d__[i__4].r;
	z__1.r = d__[i__3].r + z__2.r, z__1.i = d__[i__3].i + z__2.i;
	c__[i__2].r = z__1.r, c__[i__2].i = z__1.i;
	i__2 = kp1;
	i__3 = kp1;
	i__4 = k;
	z__2.r = t.r * e[i__4].r - t.i * e[i__4].i, z__2.i = t.r * e[i__4].i 
		+ t.i * e[i__4].r;
	z__1.r = e[i__3].r + z__2.r, z__1.i = e[i__3].i + z__2.i;
	d__[i__2].r = z__1.r, d__[i__2].i = z__1.i;
	i__2 = kp1;
	e[i__2].r = 0., e[i__2].i = 0.;
	i__2 = kp1;
	i__3 = kp1;
	i__4 = k;
	z__2.r = t.r * b[i__4].r - t.i * b[i__4].i, z__2.i = t.r * b[i__4].i 
		+ t.i * b[i__4].r;
	z__1.r = b[i__3].r + z__2.r, z__1.i = b[i__3].i + z__2.i;
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L30: */
    }
L40:
    i__1 = *n;
    if ((d__1 = c__[i__1].r, abs(d__1)) + (d__2 = d_imag(&c__[*n]), abs(d__2))
	     != 0.) {
	goto L50;
    }
    *info = *n;
    goto L90;
L50:

/*           BACK SOLVE */

    nm2 = *n - 2;
    i__1 = *n;
    z_div(&z__1, &b[*n], &c__[*n]);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    if (*n == 1) {
	goto L80;
    }
    i__1 = nm1;
    i__2 = nm1;
    i__3 = nm1;
    i__4 = *n;
    z__3.r = d__[i__3].r * b[i__4].r - d__[i__3].i * b[i__4].i, z__3.i = d__[
	    i__3].r * b[i__4].i + d__[i__3].i * b[i__4].r;
    z__2.r = b[i__2].r - z__3.r, z__2.i = b[i__2].i - z__3.i;
    z_div(&z__1, &z__2, &c__[nm1]);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    if (nm2 < 1) {
	goto L70;
    }
    i__1 = nm2;
    for (kb = 1; kb <= i__1; ++kb) {
	k = nm2 - kb + 1;
	i__2 = k;
	i__3 = k;
	i__4 = k;
	i__5 = k + 1;
	z__4.r = d__[i__4].r * b[i__5].r - d__[i__4].i * b[i__5].i, z__4.i = 
		d__[i__4].r * b[i__5].i + d__[i__4].i * b[i__5].r;
	z__3.r = b[i__3].r - z__4.r, z__3.i = b[i__3].i - z__4.i;
	i__6 = k;
	i__7 = k + 2;
	z__5.r = e[i__6].r * b[i__7].r - e[i__6].i * b[i__7].i, z__5.i = e[
		i__6].r * b[i__7].i + e[i__6].i * b[i__7].r;
	z__2.r = z__3.r - z__5.r, z__2.i = z__3.i - z__5.i;
	z_div(&z__1, &z__2, &c__[k]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L60: */
    }
L70:
L80:
L90:
L100:

    return 0;
} /* cgtsl_ */

/* Subroutine */ int cptsl_(integer *n, doublecomplex *d__, doublecomplex *e, 
	doublecomplex *b)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer k;
    static doublecomplex t1, t2;
    static integer ke, kf, kp1, nm1, kbm1, nm1d2;


/*     CPTSL GIVEN A POSITIVE DEFINITE TRIDIAGONAL MATRIX AND A RIGHT */
/*     HAND SIDE WILL FIND THE SOLUTION. */

/*     ON ENTRY */

/*        N        INTEGER */
/*                 IS THE ORDER OF THE TRIDIAGONAL MATRIX. */

/*        D        COMPLEX(N) */
/*                 IS THE DIAGONAL OF THE TRIDIAGONAL MATRIX. */
/*                 ON OUTPUT D IS DESTROYED. */

/*        E        COMPLEX(N) */
/*                 IS THE OFFDIAGONAL OF THE TRIDIAGONAL MATRIX. */
/*                 E(1) THROUGH E(N-1) SHOULD CONTAIN THE */
/*                 OFFDIAGONAL. */

/*        B        COMPLEX(N) */
/*                 IS THE RIGHT HAND SIDE VECTOR. */

/*     ON RETURN */

/*        B        CONTAINS THE SOULTION. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     JACK DONGARRA, ARGONNE NATIONAL LABORATORY. */

/*     NO EXTERNALS */
/*     FORTRAN CONJG,MOD */

/*     INTERNAL VARIABLES */


/*     CHECK FOR 1 X 1 CASE */

    /* Parameter adjustments */
    --b;
    --e;
    --d__;

    /* Function Body */
    if (*n != 1) {
	goto L10;
    }
    z_div(&z__1, &b[1], &d__[1]);
    b[1].r = z__1.r, b[1].i = z__1.i;
    goto L70;
L10:
    nm1 = *n - 1;
    nm1d2 = nm1 / 2;
    if (*n == 2) {
	goto L30;
    }
    kbm1 = *n - 1;

/*           ZERO TOP HALF OF SUBDIAGONAL AND BOTTOM HALF OF */
/*           SUPERDIAGONAL */

    i__1 = nm1d2;
    for (k = 1; k <= i__1; ++k) {
	d_cnjg(&z__2, &e[k]);
	z_div(&z__1, &z__2, &d__[k]);
	t1.r = z__1.r, t1.i = z__1.i;
	i__2 = k + 1;
	i__3 = k + 1;
	i__4 = k;
	z__2.r = t1.r * e[i__4].r - t1.i * e[i__4].i, z__2.i = t1.r * e[i__4]
		.i + t1.i * e[i__4].r;
	z__1.r = d__[i__3].r - z__2.r, z__1.i = d__[i__3].i - z__2.i;
	d__[i__2].r = z__1.r, d__[i__2].i = z__1.i;
	i__2 = k + 1;
	i__3 = k + 1;
	i__4 = k;
	z__2.r = t1.r * b[i__4].r - t1.i * b[i__4].i, z__2.i = t1.r * b[i__4]
		.i + t1.i * b[i__4].r;
	z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	z_div(&z__1, &e[kbm1], &d__[kbm1 + 1]);
	t2.r = z__1.r, t2.i = z__1.i;
	i__2 = kbm1;
	i__3 = kbm1;
	d_cnjg(&z__3, &e[kbm1]);
	z__2.r = t2.r * z__3.r - t2.i * z__3.i, z__2.i = t2.r * z__3.i + t2.i 
		* z__3.r;
	z__1.r = d__[i__3].r - z__2.r, z__1.i = d__[i__3].i - z__2.i;
	d__[i__2].r = z__1.r, d__[i__2].i = z__1.i;
	i__2 = kbm1;
	i__3 = kbm1;
	i__4 = kbm1 + 1;
	z__2.r = t2.r * b[i__4].r - t2.i * b[i__4].i, z__2.i = t2.r * b[i__4]
		.i + t2.i * b[i__4].r;
	z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	--kbm1;
/* L20: */
    }
L30:
    kp1 = nm1d2 + 1;

/*        CLEAN UP FOR POSSIBLE 2 X 2 BLOCK AT CENTER */

    if (*n % 2 != 0) {
	goto L40;
    }
    d_cnjg(&z__2, &e[kp1]);
    z_div(&z__1, &z__2, &d__[kp1]);
    t1.r = z__1.r, t1.i = z__1.i;
    i__1 = kp1 + 1;
    i__2 = kp1 + 1;
    i__3 = kp1;
    z__2.r = t1.r * e[i__3].r - t1.i * e[i__3].i, z__2.i = t1.r * e[i__3].i + 
	    t1.i * e[i__3].r;
    z__1.r = d__[i__2].r - z__2.r, z__1.i = d__[i__2].i - z__2.i;
    d__[i__1].r = z__1.r, d__[i__1].i = z__1.i;
    i__1 = kp1 + 1;
    i__2 = kp1 + 1;
    i__3 = kp1;
    z__2.r = t1.r * b[i__3].r - t1.i * b[i__3].i, z__2.i = t1.r * b[i__3].i + 
	    t1.i * b[i__3].r;
    z__1.r = b[i__2].r - z__2.r, z__1.i = b[i__2].i - z__2.i;
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    ++kp1;
L40:

/*        BACK SOLVE STARTING AT THE CENTER, GOING TOWARDS THE TOP */
/*        AND BOTTOM */

    i__1 = kp1;
    z_div(&z__1, &b[kp1], &d__[kp1]);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    if (*n == 2) {
	goto L60;
    }
    k = kp1 - 1;
    ke = kp1 + nm1d2 - 1;
    i__1 = ke;
    for (kf = kp1; kf <= i__1; ++kf) {
	i__2 = k;
	i__3 = k;
	i__4 = k;
	i__5 = k + 1;
	z__3.r = e[i__4].r * b[i__5].r - e[i__4].i * b[i__5].i, z__3.i = e[
		i__4].r * b[i__5].i + e[i__4].i * b[i__5].r;
	z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
	z_div(&z__1, &z__2, &d__[k]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	i__2 = kf + 1;
	i__3 = kf + 1;
	d_cnjg(&z__4, &e[kf]);
	i__4 = kf;
	z__3.r = z__4.r * b[i__4].r - z__4.i * b[i__4].i, z__3.i = z__4.r * b[
		i__4].i + z__4.i * b[i__4].r;
	z__2.r = b[i__3].r - z__3.r, z__2.i = b[i__3].i - z__3.i;
	z_div(&z__1, &z__2, &d__[kf + 1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	--k;
/* L50: */
    }
L60:
    if (*n % 2 == 0) {
	z__3.r = e[1].r * b[2].r - e[1].i * b[2].i, z__3.i = e[1].r * b[2].i 
		+ e[1].i * b[2].r;
	z__2.r = b[1].r - z__3.r, z__2.i = b[1].i - z__3.i;
	z_div(&z__1, &z__2, &d__[1]);
	b[1].r = z__1.r, b[1].i = z__1.i;
    }
L70:
    return 0;
} /* cptsl_ */

/* Subroutine */ int cchdc_(doublecomplex *a, integer *lda, integer *p, 
	doublecomplex *work, integer *jpvt, integer *job, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double sqrt(doublereal);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, k, l, kb, jp, pl, jt, pu, km1, kp1, plp1;
    static logical negk;
    static integer maxl;
    static doublecomplex temp;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static logical swapk;
    static doublereal maxdia;


/*     CCHDC COMPUTES THE CHOLESKY DECOMPOSITION OF A POSITIVE DEFINITE */
/*     MATRIX.  A PIVOTING OPTION ALLOWS THE USER TO ESTIMATE THE */
/*     CONDITION OF A POSITIVE DEFINITE MATRIX OR DETERMINE THE RANK */
/*     OF A POSITIVE SEMIDEFINITE MATRIX. */

/*     ON ENTRY */

/*         A      COMPLEX(LDA,P). */
/*                A CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO */
/*                BE COMPUTED.  ONLT THE UPPER HALF OF A NEED BE STORED. */
/*                THE LOWER PART OF THE ARRAY A IS NOT REFERENCED. */

/*         LDA    INTEGER. */
/*                LDA IS THE LEADING DIMENSION OF THE ARRAY A. */

/*         P      INTEGER. */
/*                P IS THE ORDER OF THE MATRIX. */

/*         WORK   COMPLEX. */
/*                WORK IS A WORK ARRAY. */

/*         JPVT   INTEGER(P). */
/*                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION */
/*                OF THE PIVOT ELEMENTS, IF PIVOTING HAS BEEN REQUESTED. */
/*                EACH DIAGONAL ELEMENT A(K,K) */
/*                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE */
/*                VALUE OF JPVT(K). */

/*                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL */
/*                                      ELEMENT. */

/*                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE ELEMENT. */

/*                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL ELEMENT. */

/*                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL ELEMENTS */
/*                ARE MOVED BY SYMMETRIC ROW AND COLUMN INTERCHANGES TO */
/*                THE BEGINNING OF THE ARRAY A AND FINAL */
/*                ELEMENTS TO THE END.  BOTH INITIAL AND FINAL ELEMENTS */
/*                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY */
/*                FREE ELEMENTS ARE MOVED.  AT THE K-TH STAGE OF THE */
/*                REDUCTION, IF A(K,K) IS OCCUPIED BY A FREE ELEMENT */
/*                IT IS INTERCHANGED WITH THE LARGEST FREE ELEMENT */
/*                A(L,L) WITH L .GE. K.  JPVT IS NOT REFERENCED IF */
/*                JOB .EQ. 0. */

/*        JOB     INTEGER. */
/*                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING. */
/*                IF JOB .EQ. 0, NO PIVOTING IS DONE. */
/*                IF JOB .NE. 0, PIVOTING IS DONE. */

/*     ON RETURN */

/*         A      A CONTAINS IN ITS UPPER HALF THE CHOLESKY FACTOR */
/*                OF THE MATRIX A AS IT HAS BEEN PERMUTED BY PIVOTING. */

/*         JPVT   JPVT(J) CONTAINS THE INDEX OF THE DIAGONAL ELEMENT */
/*                OF A THAT WAS MOVED INTO THE J-TH POSITION, */
/*                PROVIDED PIVOTING WAS REQUESTED. */

/*         INFO   CONTAINS THE INDEX OF THE LAST POSITIVE DIAGONAL */
/*                ELEMENT OF THE CHOLESKY FACTOR. */

/*     FOR POSITIVE DEFINITE MATRICES INFO = P IS THE NORMAL RETURN. */
/*     FOR PIVOTING WITH POSITIVE SEMIDEFINITE MATRICES INFO WILL */
/*     IN GENERAL BE LESS THAN P.  HOWEVER, INFO MAY BE GREATER THAN */
/*     THE RANK OF A, SINCE ROUNDING ERROR CAN CAUSE AN OTHERWISE ZERO */
/*     ELEMENT TO BE POSITIVE. INDEFINITE SYSTEMS WILL ALWAYS CAUSE */
/*     INFO TO BE LESS THAN P. */

/*     LINPACK. THIS VERSION DATED 03/19/79 . */
/*     J.J. DONGARRA AND G.W. STEWART, ARGONNE NATIONAL LABORATORY AND */
/*     UNIVERSITY OF MARYLAND. */


/*     BLAS CAXPY,CSWAP */
/*     FORTRAN SQRT,DBLE,CONJG */

/*     INTERNAL VARIABLES */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    --jpvt;

    /* Function Body */
    pl = 1;
    pu = 0;
    *info = *p;
    if (*job == 0) {
	goto L160;
    }

/*        PIVOTING HAS BEEN REQUESTED. REARRANGE THE */
/*        THE ELEMENTS ACCORDING TO JPVT. */

    i__1 = *p;
    for (k = 1; k <= i__1; ++k) {
	swapk = jpvt[k] > 0;
	negk = jpvt[k] < 0;
	jpvt[k] = k;
	if (negk) {
	    jpvt[k] = -jpvt[k];
	}
	if (! swapk) {
	    goto L60;
	}
	if (k == pl) {
	    goto L50;
	}
	i__2 = pl - 1;
	cswap_(&i__2, &a[k * a_dim1 + 1], &c__1, &a[pl * a_dim1 + 1], &c__1);
	i__2 = k + k * a_dim1;
	temp.r = a[i__2].r, temp.i = a[i__2].i;
	i__2 = k + k * a_dim1;
	i__3 = pl + pl * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = pl + pl * a_dim1;
	a[i__2].r = temp.r, a[i__2].i = temp.i;
	i__2 = pl + k * a_dim1;
	d_cnjg(&z__1, &a[pl + k * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	plp1 = pl + 1;
	if (*p < plp1) {
	    goto L40;
	}
	i__2 = *p;
	for (j = plp1; j <= i__2; ++j) {
	    if (j >= k) {
		goto L10;
	    }
	    d_cnjg(&z__1, &a[pl + j * a_dim1]);
	    temp.r = z__1.r, temp.i = z__1.i;
	    i__3 = pl + j * a_dim1;
	    d_cnjg(&z__1, &a[j + k * a_dim1]);
	    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
	    i__3 = j + k * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
	    goto L20;
L10:
	    if (j == k) {
		goto L20;
	    }
	    i__3 = k + j * a_dim1;
	    temp.r = a[i__3].r, temp.i = a[i__3].i;
	    i__3 = k + j * a_dim1;
	    i__4 = pl + j * a_dim1;
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
	    i__3 = pl + j * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
L20:
/* L30: */
	    ;
	}
L40:
	jpvt[k] = jpvt[pl];
	jpvt[pl] = k;
L50:
	++pl;
L60:
/* L70: */
	;
    }
    pu = *p;
    if (*p < pl) {
	goto L150;
    }
    i__1 = *p;
    for (kb = pl; kb <= i__1; ++kb) {
	k = *p - kb + pl;
	if (jpvt[k] >= 0) {
	    goto L130;
	}
	jpvt[k] = -jpvt[k];
	if (pu == k) {
	    goto L120;
	}
	i__2 = k - 1;
	cswap_(&i__2, &a[k * a_dim1 + 1], &c__1, &a[pu * a_dim1 + 1], &c__1);
	i__2 = k + k * a_dim1;
	temp.r = a[i__2].r, temp.i = a[i__2].i;
	i__2 = k + k * a_dim1;
	i__3 = pu + pu * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = pu + pu * a_dim1;
	a[i__2].r = temp.r, a[i__2].i = temp.i;
	i__2 = k + pu * a_dim1;
	d_cnjg(&z__1, &a[k + pu * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	kp1 = k + 1;
	if (*p < kp1) {
	    goto L110;
	}
	i__2 = *p;
	for (j = kp1; j <= i__2; ++j) {
	    if (j >= pu) {
		goto L80;
	    }
	    d_cnjg(&z__1, &a[k + j * a_dim1]);
	    temp.r = z__1.r, temp.i = z__1.i;
	    i__3 = k + j * a_dim1;
	    d_cnjg(&z__1, &a[j + pu * a_dim1]);
	    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
	    i__3 = j + pu * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
	    goto L90;
L80:
	    if (j == pu) {
		goto L90;
	    }
	    i__3 = k + j * a_dim1;
	    temp.r = a[i__3].r, temp.i = a[i__3].i;
	    i__3 = k + j * a_dim1;
	    i__4 = pu + j * a_dim1;
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
	    i__3 = pu + j * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
L90:
/* L100: */
	    ;
	}
L110:
	jt = jpvt[k];
	jpvt[k] = jpvt[pu];
	jpvt[pu] = jt;
L120:
	--pu;
L130:
/* L140: */
	;
    }
L150:
L160:
    i__1 = *p;
    for (k = 1; k <= i__1; ++k) {

/*        REDUCTION LOOP. */

	i__2 = k + k * a_dim1;
	maxdia = a[i__2].r;
	kp1 = k + 1;
	maxl = k;

/*        DETERMINE THE PIVOT ELEMENT. */

	if (k < pl || k >= pu) {
	    goto L190;
	}
	i__2 = pu;
	for (l = kp1; l <= i__2; ++l) {
	    i__3 = l + l * a_dim1;
	    if (a[i__3].r <= maxdia) {
		goto L170;
	    }
	    i__3 = l + l * a_dim1;
	    maxdia = a[i__3].r;
	    maxl = l;
L170:
/* L180: */
	    ;
	}
L190:

/*        QUIT IF THE PIVOT ELEMENT IS NOT POSITIVE. */

	if (maxdia > 0.) {
	    goto L200;
	}
	*info = k - 1;
/*     ......EXIT */
	goto L280;
L200:
	if (k == maxl) {
	    goto L210;
	}

/*           START THE PIVOTING AND UPDATE JPVT. */

	km1 = k - 1;
	cswap_(&km1, &a[k * a_dim1 + 1], &c__1, &a[maxl * a_dim1 + 1], &c__1);
	i__2 = maxl + maxl * a_dim1;
	i__3 = k + k * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = k + k * a_dim1;
	z__1.r = maxdia, z__1.i = 0.;
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	jp = jpvt[maxl];
	jpvt[maxl] = jpvt[k];
	jpvt[k] = jp;
	i__2 = k + maxl * a_dim1;
	d_cnjg(&z__1, &a[k + maxl * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
L210:

/*        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS. */

	i__2 = k;
	d__1 = sqrt((doublereal) a[k + k * a_dim1].r);
	z__1.r = d__1, z__1.i = 0.;
	work[i__2].r = z__1.r, work[i__2].i = z__1.i;
	i__2 = k + k * a_dim1;
	i__3 = k;
	a[i__2].r = work[i__3].r, a[i__2].i = work[i__3].i;
	if (*p < kp1) {
	    goto L260;
	}
	i__2 = *p;
	for (j = kp1; j <= i__2; ++j) {
	    if (k == maxl) {
		goto L240;
	    }
	    if (j >= maxl) {
		goto L220;
	    }
	    d_cnjg(&z__1, &a[k + j * a_dim1]);
	    temp.r = z__1.r, temp.i = z__1.i;
	    i__3 = k + j * a_dim1;
	    d_cnjg(&z__1, &a[j + maxl * a_dim1]);
	    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
	    i__3 = j + maxl * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
	    goto L230;
L220:
	    if (j == maxl) {
		goto L230;
	    }
	    i__3 = k + j * a_dim1;
	    temp.r = a[i__3].r, temp.i = a[i__3].i;
	    i__3 = k + j * a_dim1;
	    i__4 = maxl + j * a_dim1;
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
	    i__3 = maxl + j * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
L230:
L240:
	    i__3 = k + j * a_dim1;
	    z_div(&z__1, &a[k + j * a_dim1], &work[k]);
	    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
	    i__3 = j;
	    d_cnjg(&z__1, &a[k + j * a_dim1]);
	    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
	    i__3 = k + j * a_dim1;
	    z__1.r = -a[i__3].r, z__1.i = -a[i__3].i;
	    temp.r = z__1.r, temp.i = z__1.i;
	    i__3 = j - k;
	    caxpy_(&i__3, &temp, &work[kp1], &c__1, &a[kp1 + j * a_dim1], &
		    c__1);
/* L250: */
	}
L260:
/* L270: */
	;
    }
L280:
    return 0;
} /* cchdc_ */

/* Subroutine */ int cchud_(doublecomplex *r__, integer *ldr, integer *p, 
	doublecomplex *x, doublecomplex *z__, integer *ldz, integer *nz, 
	doublecomplex *y, doublereal *rho, doublereal *c__, doublecomplex *s)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double z_abs(doublecomplex *), sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublecomplex t, xj;
    static integer jm1;
    static doublecomplex zeta;
    static doublereal scale, azeta;
    extern /* Subroutine */ int crotg_(doublecomplex *, doublecomplex *, 
	    doublereal *, doublecomplex *);


/*     CCHUD UPDATES AN AUGMENTED CHOLESKY DECOMPOSITION OF THE */
/*     TRIANGULAR PART OF AN AUGMENTED QR DECOMPOSITION.  SPECIFICALLY, */
/*     GIVEN AN UPPER TRIANGULAR MATRIX R OF ORDER P, A ROW VECTOR */
/*     X, A COLUMN VECTOR Z, AND A SCALAR Y, CCHUD DETERMINES A */
/*     UNTIARY MATRIX U AND A SCALAR ZETA SUCH THAT */


/*                              (R  Z)     (RR   ZZ ) */
/*                         U  * (    )  =  (        ) , */
/*                              (X  Y)     ( 0  ZETA) */

/*     WHERE RR IS UPPER TRIANGULAR.  IF R AND Z HAVE BEEN */
/*     OBTAINED FROM THE FACTORIZATION OF A LEAST SQUARES */
/*     PROBLEM, THEN RR AND ZZ ARE THE FACTORS CORRESPONDING TO */
/*     THE PROBLEM WITH THE OBSERVATION (X,Y) APPENDED.  IN THIS */
/*     CASE, IF RHO IS THE NORM OF THE RESIDUAL VECTOR, THEN THE */
/*     NORM OF THE RESIDUAL VECTOR OF THE UPDATED PROBLEM IS */
/*     SQRT(RHO**2 + ZETA**2).  CCHUD WILL SIMULTANEOUSLY UPDATE */
/*     SEVERAL TRIPLETS (Z,Y,RHO). */
/*     FOR A LESS TERSE DESCRIPTION OF WHAT CCHUD DOES AND HOW */
/*     IT MAY BE APPLIED SEE THE LINPACK GUIDE. */

/*     THE MATRIX U IS DETERMINED AS THE PRODUCT U(P)*...*U(1), */
/*     WHERE U(I) IS A ROTATION IN THE (I,P+1) PLANE OF THE */
/*     FORM */

/*                       (     C(I)      S(I) ) */
/*                       (                    ) . */
/*                       ( -CONJG(S(I))  C(I) ) */

/*     THE ROTATIONS ARE CHOSEN SO THAT C(I) IS REAL. */

/*     ON ENTRY */

/*         R      COMPLEX(LDR,P), WHERE LDR .GE. P. */
/*                R CONTAINS THE UPPER TRIANGULAR MATRIX */
/*                THAT IS TO BE UPDATED.  THE PART OF R */
/*                BELOW THE DIAGONAL IS NOT REFERENCED. */

/*         LDR    INTEGER. */
/*                LDR IS THE LEADING DIMENSION OF THE ARRAY R. */

/*         P      INTEGER. */
/*                P IS THE ORDER OF THE MATRIX R. */

/*         X      COMPLEX(P). */
/*                X CONTAINS THE ROW TO BE ADDED TO R.  X IS */
/*                NOT ALTERED BY CCHUD. */

/*         Z      COMPLEX(LDZ,NZ), WHERE LDZ .GE. P. */
/*                Z IS AN ARRAY CONTAINING NZ P-VECTORS TO */
/*                BE UPDATED WITH R. */

/*         LDZ    INTEGER. */
/*                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z. */

/*         NZ     INTEGER. */
/*                NZ IS THE NUMBER OF VECTORS TO BE UPDATED */
/*                NZ MAY BE ZERO, IN WHICH CASE Z, Y, AND RHO */
/*                ARE NOT REFERENCED. */

/*         Y      COMPLEX(NZ). */
/*                Y CONTAINS THE SCALARS FOR UPDATING THE VECTORS */
/*                Z.  Y IS NOT ALTERED BY CCHUD. */

/*         RHO    REAL(NZ). */
/*                RHO CONTAINS THE NORMS OF THE RESIDUAL */
/*                VECTORS THAT ARE TO BE UPDATED.  IF RHO(J) */
/*                IS NEGATIVE, IT IS LEFT UNALTERED. */

/*     ON RETURN */

/*         RC */
/*         RHO    CONTAIN THE UPDATED QUANTITIES. */
/*         Z */

/*         C      REAL(P). */
/*                C CONTAINS THE COSINES OF THE TRANSFORMING */
/*                ROTATIONS. */

/*         S      COMPLEX(P). */
/*                S CONTAINS THE SINES OF THE TRANSFORMING */
/*                ROTATIONS. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     CCHUD USES THE FOLLOWING FUNCTIONS AND SUBROUTINES. */

/*     EXTENDED BLAS CROTG */
/*     FORTRAN CONJG,SQRT */


/*     UPDATE R. */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --x;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --y;
    --rho;
    --c__;
    --s;

    /* Function Body */
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	xj.r = x[i__2].r, xj.i = x[i__2].i;

/*        APPLY THE PREVIOUS ROTATIONS. */

	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__ + j * r_dim1;
	    z__2.r = c__[i__3] * r__[i__4].r, z__2.i = c__[i__3] * r__[i__4]
		    .i;
	    i__5 = i__;
	    z__3.r = s[i__5].r * xj.r - s[i__5].i * xj.i, z__3.i = s[i__5].r *
		     xj.i + s[i__5].i * xj.r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = i__;
	    z__2.r = c__[i__3] * xj.r, z__2.i = c__[i__3] * xj.i;
	    d_cnjg(&z__4, &s[i__]);
	    i__4 = i__ + j * r_dim1;
	    z__3.r = z__4.r * r__[i__4].r - z__4.i * r__[i__4].i, z__3.i = 
		    z__4.r * r__[i__4].i + z__4.i * r__[i__4].r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    xj.r = z__1.r, xj.i = z__1.i;
	    i__3 = i__ + j * r_dim1;
	    r__[i__3].r = t.r, r__[i__3].i = t.i;
/* L10: */
	}
L20:

/*        COMPUTE THE NEXT ROTATION. */

	crotg_(&r__[j + j * r_dim1], &xj, &c__[j], &s[j]);
/* L30: */
    }

/*     IF REQUIRED, UPDATE Z AND RHO. */

    if (*nz < 1) {
	goto L70;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	zeta.r = y[i__2].r, zeta.i = y[i__2].i;
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__ + j * z_dim1;
	    z__2.r = c__[i__3] * z__[i__4].r, z__2.i = c__[i__3] * z__[i__4]
		    .i;
	    i__5 = i__;
	    z__3.r = s[i__5].r * zeta.r - s[i__5].i * zeta.i, z__3.i = s[i__5]
		    .r * zeta.i + s[i__5].i * zeta.r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = i__;
	    z__2.r = c__[i__3] * zeta.r, z__2.i = c__[i__3] * zeta.i;
	    d_cnjg(&z__4, &s[i__]);
	    i__4 = i__ + j * z_dim1;
	    z__3.r = z__4.r * z__[i__4].r - z__4.i * z__[i__4].i, z__3.i = 
		    z__4.r * z__[i__4].i + z__4.i * z__[i__4].r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    zeta.r = z__1.r, zeta.i = z__1.i;
	    i__3 = i__ + j * z_dim1;
	    z__[i__3].r = t.r, z__[i__3].i = t.i;
/* L40: */
	}
	azeta = z_abs(&zeta);
	if (azeta == 0. || rho[j] < 0.) {
	    goto L50;
	}
	scale = azeta + rho[j];
/* Computing 2nd power */
	d__1 = azeta / scale;
/* Computing 2nd power */
	d__2 = rho[j] / scale;
	rho[j] = scale * sqrt(d__1 * d__1 + d__2 * d__2);
L50:
/* L60: */
	;
    }
L70:
    return 0;
} /* cchud_ */

/* Subroutine */ int cchdd_(doublecomplex *r__, integer *ldr, integer *p, 
	doublecomplex *x, doublecomplex *z__, integer *ldz, integer *nz, 
	doublecomplex *y, doublereal *rho, doublereal *c__, doublecomplex *s, 
	integer *info)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);
    double sqrt(doublereal), z_abs(doublecomplex *), d_imag(doublecomplex *);

    /* Local variables */
    static doublereal a;
    static doublecomplex b;
    static integer i__, j;
    static doublecomplex t;
    static integer ii;
    static doublecomplex xx, zeta;
    static doublereal norm, alpha;
    static real scale;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal azeta;
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);


/*     CCHDD DOWNDATES AN AUGMENTED CHOLESKY DECOMPOSITION OR THE */
/*     TRIANGULAR FACTOR OF AN AUGMENTED QR DECOMPOSITION. */
/*     SPECIFICALLY, GIVEN AN UPPER TRIANGULAR MATRIX R OF ORDER P,  A */
/*     ROW VECTOR X, A COLUMN VECTOR Z, AND A SCALAR Y, CCHDD */
/*     DETERMINEDS A UNITARY MATRIX U AND A SCALAR ZETA SUCH THAT */

/*                        (R   Z )     (RR  ZZ) */
/*                    U * (      )  =  (      ) , */
/*                        (0 ZETA)     ( X   Y) */

/*     WHERE RR IS UPPER TRIANGULAR.  IF R AND Z HAVE BEEN OBTAINED */
/*     FROM THE FACTORIZATION OF A LEAST SQUARES PROBLEM, THEN */
/*     RR AND ZZ ARE THE FACTORS CORRESPONDING TO THE PROBLEM */
/*     WITH THE OBSERVATION (X,Y) REMOVED.  IN THIS CASE, IF RHO */
/*     IS THE NORM OF THE RESIDUAL VECTOR, THEN THE NORM OF */
/*     THE RESIDUAL VECTOR OF THE DOWNDATED PROBLEM IS */
/*     SQRT(RHO**2 - ZETA**2). CCHDD WILL SIMULTANEOUSLY DOWNDATE */
/*     SEVERAL TRIPLETS (Z,Y,RHO) ALONG WITH R. */
/*     FOR A LESS TERSE DESCRIPTION OF WHAT CCHDD DOES AND HOW */
/*     IT MAY BE APPLIED, SEE THE LINPACK GUIDE. */

/*     THE MATRIX U IS DETERMINED AS THE PRODUCT U(1)*...*U(P) */
/*     WHERE U(I) IS A ROTATION IN THE (P+1,I)-PLANE OF THE */
/*     FORM */

/*                       ( C(I)  -CONJG(S(I)) ) */
/*                       (                    ) . */
/*                       ( S(I)       C(I)    ) */

/*     THE ROTATIONS ARE CHOSEN SO THAT C(I) IS REAL. */

/*     THE USER IS WARNED THAT A GIVEN DOWNDATING PROBLEM MAY */
/*     BE IMPOSSIBLE TO ACCOMPLISH OR MAY PRODUCE */
/*     INACCURATE RESULTS.  FOR EXAMPLE, THIS CAN HAPPEN */
/*     IF X IS NEAR A VECTOR WHOSE REMOVAL WILL REDUCE THE */
/*     RANK OF R.  BEWARE. */

/*     ON ENTRY */

/*         R      COMPLEX(LDR,P), WHERE LDR .GE. P. */
/*                R CONTAINS THE UPPER TRIANGULAR MATRIX */
/*                THAT IS TO BE DOWNDATED.  THE PART OF  R */
/*                BELOW THE DIAGONAL IS NOT REFERENCED. */

/*         LDR    INTEGER. */
/*                LDR IS THE LEADING DIMENSION FO THE ARRAY R. */

/*         P      INTEGER. */
/*                P IS THE ORDER OF THE MATRIX R. */

/*         X      COMPLEX(P). */
/*                X CONTAINS THE ROW VECTOR THAT IS TO */
/*                BE REMOVED FROM R.  X IS NOT ALTERED BY CCHDD. */

/*         Z      COMPLEX(LDZ,NZ), WHERE LDZ .GE. P. */
/*                Z IS AN ARRAY OF NZ P-VECTORS WHICH */
/*                ARE TO BE DOWNDATED ALONG WITH R. */

/*         LDZ    INTEGER. */
/*                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z. */

/*         NZ     INTEGER. */
/*                NZ IS THE NUMBER OF VECTORS TO BE DOWNDATED */
/*                NZ MAY BE ZERO, IN WHICH CASE Z, Y, AND RHO */
/*                ARE NOT REFERENCED. */

/*         Y      COMPLEX(NZ). */
/*                Y CONTAINS THE SCALARS FOR THE DOWNDATING */
/*                OF THE VECTORS Z.  Y IS NOT ALTERED BY CCHDD. */

/*         RHO    REAL(NZ). */
/*                RHO CONTAINS THE NORMS OF THE RESIDUAL */
/*                VECTORS THAT ARE TO BE DOWNDATED. */

/*     ON RETURN */

/*         R */
/*         Z      CONTAIN THE DOWNDATED QUANTITIES. */
/*         RHO */

/*         C      REAL(P). */
/*                C CONTAINS THE COSINES OF THE TRANSFORMING */
/*                ROTATIONS. */

/*         S      COMPLEX(P). */
/*                S CONTAINS THE SINES OF THE TRANSFORMING */
/*                ROTATIONS. */

/*         INFO   INTEGER. */
/*                INFO IS SET AS FOLLOWS. */

/*                   INFO = 0  IF THE ENTIRE DOWNDATING */
/*                             WAS SUCCESSFUL. */

/*                   INFO =-1  IF R COULD NOT BE DOWNDATED. */
/*                             IN THIS CASE, ALL QUANTITIES */
/*                             ARE LEFT UNALTERED. */

/*                   INFO = 1  IF SOME RHO COULD NOT BE */
/*                             DOWNDATED.  THE OFFENDING RHOS ARE */
/*                             SET TO -1. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     CCHDD USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS. */

/*     FORTRAN DIMAG,CABS,CONJG */
/*     BLAS CDOTC, SCNRM2 */


/*     SOLVE THE SYSTEM CTRANS(R)*A = X, PLACING THE RESULT */
/*     IN THE ARRAY S. */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --x;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --y;
    --rho;
    --c__;
    --s;

    /* Function Body */
    *info = 0;
    d_cnjg(&z__2, &x[1]);
    d_cnjg(&z__3, &r__[r_dim1 + 1]);
    z_div(&z__1, &z__2, &z__3);
    s[1].r = z__1.r, s[1].i = z__1.i;
    if (*p < 2) {
	goto L20;
    }
    i__1 = *p;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	d_cnjg(&z__2, &x[j]);
	i__3 = j - 1;
	cdotc_(&z__3, &i__3, &r__[j * r_dim1 + 1], &c__1, &s[1], &c__1);
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
	i__2 = j;
	d_cnjg(&z__2, &r__[j + j * r_dim1]);
	z_div(&z__1, &s[j], &z__2);
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
/* L10: */
    }
L20:
    norm = scnrm2_(p, &s[1], &c__1);
    if (norm < 1.) {
	goto L30;
    }
    *info = -1;
    goto L120;
L30:
/* Computing 2nd power */
    d__1 = norm;
    alpha = sqrt(1. - d__1 * d__1);

/*        DETERMINE THE TRANSFORMATIONS. */

    i__1 = *p;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *p - ii + 1;
	scale = alpha + z_abs(&s[i__]);
	a = alpha / scale;
	i__2 = i__;
	z__1.r = s[i__2].r / scale, z__1.i = s[i__2].i / scale;
	b.r = z__1.r, b.i = z__1.i;
/* Computing 2nd power */
	d__1 = a;
/* Computing 2nd power */
	d__2 = b.r;
/* Computing 2nd power */
	d__3 = d_imag(&b);
	norm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	c__[i__] = a / norm;
	i__2 = i__;
	d_cnjg(&z__2, &b);
	z__1.r = z__2.r / norm, z__1.i = z__2.i / norm;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
	alpha = scale * norm;
/* L40: */
    }

/*        APPLY THE TRANSFORMATIONS TO R. */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	xx.r = 0., xx.i = 0.;
	i__2 = j;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = j - ii + 1;
	    i__3 = i__;
	    z__2.r = c__[i__3] * xx.r, z__2.i = c__[i__3] * xx.i;
	    i__4 = i__;
	    i__5 = i__ + j * r_dim1;
	    z__3.r = s[i__4].r * r__[i__5].r - s[i__4].i * r__[i__5].i, 
		    z__3.i = s[i__4].r * r__[i__5].i + s[i__4].i * r__[i__5]
		    .r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = i__ + j * r_dim1;
	    i__4 = i__;
	    i__5 = i__ + j * r_dim1;
	    z__2.r = c__[i__4] * r__[i__5].r, z__2.i = c__[i__4] * r__[i__5]
		    .i;
	    d_cnjg(&z__4, &s[i__]);
	    z__3.r = z__4.r * xx.r - z__4.i * xx.i, z__3.i = z__4.r * xx.i + 
		    z__4.i * xx.r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    r__[i__3].r = z__1.r, r__[i__3].i = z__1.i;
	    xx.r = t.r, xx.i = t.i;
/* L50: */
	}
/* L60: */
    }

/*        IF REQUIRED, DOWNDATE Z AND RHO. */

    if (*nz < 1) {
	goto L110;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	zeta.r = y[i__2].r, zeta.i = y[i__2].i;
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * z_dim1;
	    i__4 = i__ + j * z_dim1;
	    d_cnjg(&z__4, &s[i__]);
	    z__3.r = z__4.r * zeta.r - z__4.i * zeta.i, z__3.i = z__4.r * 
		    zeta.i + z__4.i * zeta.r;
	    z__2.r = z__[i__4].r - z__3.r, z__2.i = z__[i__4].i - z__3.i;
	    i__5 = i__;
	    z__1.r = z__2.r / c__[i__5], z__1.i = z__2.i / c__[i__5];
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = i__;
	    z__2.r = c__[i__3] * zeta.r, z__2.i = c__[i__3] * zeta.i;
	    i__4 = i__;
	    i__5 = i__ + j * z_dim1;
	    z__3.r = s[i__4].r * z__[i__5].r - s[i__4].i * z__[i__5].i, 
		    z__3.i = s[i__4].r * z__[i__5].i + s[i__4].i * z__[i__5]
		    .r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    zeta.r = z__1.r, zeta.i = z__1.i;
/* L70: */
	}
	azeta = z_abs(&zeta);
	if (azeta <= rho[j]) {
	    goto L80;
	}
	*info = 1;
	rho[j] = -1.;
	goto L90;
L80:
/* Computing 2nd power */
	d__1 = azeta / rho[j];
	rho[j] *= sqrt(1. - d__1 * d__1);
L90:
/* L100: */
	;
    }
L110:
L120:
    return 0;
} /* cchdd_ */

/* Subroutine */ int cchex_(doublecomplex *r__, integer *ldr, integer *p, 
	integer *k, integer *l, doublecomplex *z__, integer *ldz, integer *nz,
	 doublereal *c__, doublecomplex *s, integer *job)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublecomplex t;
    static integer ii, jj, il, iu, km1, lm1, kp1, lmk;
    extern /* Subroutine */ int crotg_(doublecomplex *, doublecomplex *, 
	    doublereal *, doublecomplex *);


/*     CCHEX UPDATES THE CHOLESKY FACTORIZATION */

/*                   A = CTRANS(R)*R */

/*     OF A POSITIVE DEFINITE MATRIX A OF ORDER P UNDER DIAGONAL */
/*     PERMUTATIONS OF THE FORM */

/*                   TRANS(E)*A*E */

/*     WHERE E IS A PERMUTATION MATRIX.  SPECIFICALLY, GIVEN */
/*     AN UPPER TRIANGULAR MATRIX R AND A PERMUTATION MATRIX */
/*     E (WHICH IS SPECIFIED BY K, L, AND JOB), CCHEX DETERMINES */
/*     A UNITARY MATRIX U SUCH THAT */

/*                           U*R*E = RR, */

/*     WHERE RR IS UPPER TRIANGULAR.  AT THE USERS OPTION, THE */
/*     TRANSFORMATION U WILL BE MULTIPLIED INTO THE ARRAY Z. */
/*     IF A = CTRANS(X)*X, SO THAT R IS THE TRIANGULAR PART OF THE */
/*     QR FACTORIZATION OF X, THEN RR IS THE TRIANGULAR PART OF THE */
/*     QR FACTORIZATION OF X*E, I.E. X WITH ITS COLUMNS PERMUTED. */
/*     FOR A LESS TERSE DESCRIPTION OF WHAT CCHEX DOES AND HOW */
/*     IT MAY BE APPLIED, SEE THE LINPACK GUIDE. */

/*     THE MATRIX Q IS DETERMINED AS THE PRODUCT U(L-K)*...*U(1) */
/*     OF PLANE ROTATIONS OF THE FORM */

/*                           (    C(I)       S(I) ) */
/*                           (                    ) , */
/*                           ( -CONJG(S(I))  C(I) ) */

/*     WHERE C(I) IS REAL, THE ROWS THESE ROTATIONS OPERATE ON */
/*     ARE DESCRIBED BELOW. */

/*     THERE ARE TWO TYPES OF PERMUTATIONS, WHICH ARE DETERMINED */
/*     BY THE VALUE OF JOB. */

/*     1. RIGHT CIRCULAR SHIFT (JOB = 1). */

/*         THE COLUMNS ARE REARRANGED IN THE FOLLOWING ORDER. */

/*                1,...,K-1,L,K,K+1,...,L-1,L+1,...,P. */

/*         U IS THE PRODUCT OF L-K ROTATIONS U(I), WHERE U(I) */
/*         ACTS IN THE (L-I,L-I+1)-PLANE. */

/*     2. LEFT CIRCULAR SHIFT (JOB = 2). */
/*         THE COLUMNS ARE REARRANGED IN THE FOLLOWING ORDER */

/*                1,...,K-1,K+1,K+2,...,L,K,L+1,...,P. */

/*         U IS THE PRODUCT OF L-K ROTATIONS U(I), WHERE U(I) */
/*         ACTS IN THE (K+I-1,K+I)-PLANE. */

/*     ON ENTRY */

/*         R      COMPLEX(LDR,P), WHERE LDR.GE.P. */
/*                R CONTAINS THE UPPER TRIANGULAR FACTOR */
/*                THAT IS TO BE UPDATED.  ELEMENTS OF R */
/*                BELOW THE DIAGONAL ARE NOT REFERENCED. */

/*         LDR    INTEGER. */
/*                LDR IS THE LEADING DIMENSION OF THE ARRAY R. */

/*         P      INTEGER. */
/*                P IS THE ORDER OF THE MATRIX R. */

/*         K      INTEGER. */
/*                K IS THE FIRST COLUMN TO BE PERMUTED. */

/*         L      INTEGER. */
/*                L IS THE LAST COLUMN TO BE PERMUTED. */
/*                L MUST BE STRICTLY GREATER THAN K. */

/*         Z      COMPLEX(LDZ,NZ), WHERE LDZ.GE.P. */
/*                Z IS AN ARRAY OF NZ P-VECTORS INTO WHICH THE */
/*                TRANSFORMATION U IS MULTIPLIED.  Z IS */
/*                NOT REFERENCED IF NZ = 0. */

/*         LDZ    INTEGER. */
/*                LDZ IS THE LEADING DIMENSION OF THE ARRAY Z. */

/*         NZ     INTEGER. */
/*                NZ IS THE NUMBER OF COLUMNS OF THE MATRIX Z. */

/*         JOB    INTEGER. */
/*                JOB DETERMINES THE TYPE OF PERMUTATION. */
/*                       JOB = 1  RIGHT CIRCULAR SHIFT. */
/*                       JOB = 2  LEFT CIRCULAR SHIFT. */

/*     ON RETURN */

/*         R      CONTAINS THE UPDATED FACTOR. */

/*         Z      CONTAINS THE UPDATED MATRIX Z. */

/*         C      REAL(P). */
/*                C CONTAINS THE COSINES OF THE TRANSFORMING ROTATIONS. */

/*         S      COMPLEX(P). */
/*                S CONTAINS THE SINES OF THE TRANSFORMING ROTATIONS. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     CCHEX USES THE FOLLOWING FUNCTIONS AND SUBROUTINES. */

/*     EXTENDED BLAS CROTG */
/*     FORTRAN MIN0 */


/*     INITIALIZE */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --c__;
    --s;

    /* Function Body */
    km1 = *k - 1;
    kp1 = *k + 1;
    lmk = *l - *k;
    lm1 = *l - 1;

/*     PERFORM THE APPROPRIATE TASK. */

    switch (*job) {
	case 1:  goto L10;
	case 2:  goto L130;
    }

/*     RIGHT CIRCULAR SHIFT. */

L10:

/*        REORDER THE COLUMNS. */

    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *l - i__ + 1;
	i__2 = i__;
	i__3 = ii + *l * r_dim1;
	s[i__2].r = r__[i__3].r, s[i__2].i = r__[i__3].i;
/* L20: */
    }
    i__1 = lm1;
    for (jj = *k; jj <= i__1; ++jj) {
	j = lm1 - jj + *k;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + (j + 1) * r_dim1;
	    i__4 = i__ + j * r_dim1;
	    r__[i__3].r = r__[i__4].r, r__[i__3].i = r__[i__4].i;
/* L30: */
	}
	i__2 = j + 1 + (j + 1) * r_dim1;
	r__[i__2].r = 0., r__[i__2].i = 0.;
/* L40: */
    }
    if (*k == 1) {
	goto L60;
    }
    i__1 = km1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *l - i__ + 1;
	i__2 = i__ + *k * r_dim1;
	i__3 = ii;
	r__[i__2].r = s[i__3].r, r__[i__2].i = s[i__3].i;
/* L50: */
    }
L60:

/*        CALCULATE THE ROTATIONS. */

    t.r = s[1].r, t.i = s[1].i;
    i__1 = lmk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	crotg_(&s[i__ + 1], &t, &c__[i__], &s[i__]);
	i__2 = i__ + 1;
	t.r = s[i__2].r, t.i = s[i__2].i;
/* L70: */
    }
    i__1 = *k + *k * r_dim1;
    r__[i__1].r = t.r, r__[i__1].i = t.i;
    i__1 = *p;
    for (j = kp1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = 1, i__3 = *l - j + 1;
	il = max(i__2,i__3);
	i__2 = lmk;
	for (ii = il; ii <= i__2; ++ii) {
	    i__ = *l - ii;
	    i__3 = ii;
	    i__4 = i__ + j * r_dim1;
	    z__2.r = c__[i__3] * r__[i__4].r, z__2.i = c__[i__3] * r__[i__4]
		    .i;
	    i__5 = ii;
	    i__6 = i__ + 1 + j * r_dim1;
	    z__3.r = s[i__5].r * r__[i__6].r - s[i__5].i * r__[i__6].i, 
		    z__3.i = s[i__5].r * r__[i__6].i + s[i__5].i * r__[i__6]
		    .r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = i__ + 1 + j * r_dim1;
	    i__4 = ii;
	    i__5 = i__ + 1 + j * r_dim1;
	    z__2.r = c__[i__4] * r__[i__5].r, z__2.i = c__[i__4] * r__[i__5]
		    .i;
	    d_cnjg(&z__4, &s[ii]);
	    i__6 = i__ + j * r_dim1;
	    z__3.r = z__4.r * r__[i__6].r - z__4.i * r__[i__6].i, z__3.i = 
		    z__4.r * r__[i__6].i + z__4.i * r__[i__6].r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    r__[i__3].r = z__1.r, r__[i__3].i = z__1.i;
	    i__3 = i__ + j * r_dim1;
	    r__[i__3].r = t.r, r__[i__3].i = t.i;
/* L80: */
	}
/* L90: */
    }

/*        IF REQUIRED, APPLY THE TRANSFORMATIONS TO Z. */

    if (*nz < 1) {
	goto L120;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = lmk;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = *l - ii;
	    i__3 = ii;
	    i__4 = i__ + j * z_dim1;
	    z__2.r = c__[i__3] * z__[i__4].r, z__2.i = c__[i__3] * z__[i__4]
		    .i;
	    i__5 = ii;
	    i__6 = i__ + 1 + j * z_dim1;
	    z__3.r = s[i__5].r * z__[i__6].r - s[i__5].i * z__[i__6].i, 
		    z__3.i = s[i__5].r * z__[i__6].i + s[i__5].i * z__[i__6]
		    .r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = i__ + 1 + j * z_dim1;
	    i__4 = ii;
	    i__5 = i__ + 1 + j * z_dim1;
	    z__2.r = c__[i__4] * z__[i__5].r, z__2.i = c__[i__4] * z__[i__5]
		    .i;
	    d_cnjg(&z__4, &s[ii]);
	    i__6 = i__ + j * z_dim1;
	    z__3.r = z__4.r * z__[i__6].r - z__4.i * z__[i__6].i, z__3.i = 
		    z__4.r * z__[i__6].i + z__4.i * z__[i__6].r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = i__ + j * z_dim1;
	    z__[i__3].r = t.r, z__[i__3].i = t.i;
/* L100: */
	}
/* L110: */
    }
L120:
    goto L260;

/*     LEFT CIRCULAR SHIFT */

L130:

/*        REORDER THE COLUMNS */

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = lmk + i__;
	i__2 = ii;
	i__3 = i__ + *k * r_dim1;
	s[i__2].r = r__[i__3].r, s[i__2].i = r__[i__3].i;
/* L140: */
    }
    i__1 = lm1;
    for (j = *k; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * r_dim1;
	    i__4 = i__ + (j + 1) * r_dim1;
	    r__[i__3].r = r__[i__4].r, r__[i__3].i = r__[i__4].i;
/* L150: */
	}
	jj = j - km1;
	i__2 = jj;
	i__3 = j + 1 + (j + 1) * r_dim1;
	s[i__2].r = r__[i__3].r, s[i__2].i = r__[i__3].i;
/* L160: */
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = lmk + i__;
	i__2 = i__ + *l * r_dim1;
	i__3 = ii;
	r__[i__2].r = s[i__3].r, r__[i__2].i = s[i__3].i;
/* L170: */
    }
    i__1 = *l;
    for (i__ = kp1; i__ <= i__1; ++i__) {
	i__2 = i__ + *l * r_dim1;
	r__[i__2].r = 0., r__[i__2].i = 0.;
/* L180: */
    }

/*        REDUCTION LOOP. */

    i__1 = *p;
    for (j = *k; j <= i__1; ++j) {
	if (j == *k) {
	    goto L200;
	}

/*              APPLY THE ROTATIONS. */

/* Computing MIN */
	i__2 = j - 1, i__3 = *l - 1;
	iu = min(i__2,i__3);
	i__2 = iu;
	for (i__ = *k; i__ <= i__2; ++i__) {
	    ii = i__ - *k + 1;
	    i__3 = ii;
	    i__4 = i__ + j * r_dim1;
	    z__2.r = c__[i__3] * r__[i__4].r, z__2.i = c__[i__3] * r__[i__4]
		    .i;
	    i__5 = ii;
	    i__6 = i__ + 1 + j * r_dim1;
	    z__3.r = s[i__5].r * r__[i__6].r - s[i__5].i * r__[i__6].i, 
		    z__3.i = s[i__5].r * r__[i__6].i + s[i__5].i * r__[i__6]
		    .r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = i__ + 1 + j * r_dim1;
	    i__4 = ii;
	    i__5 = i__ + 1 + j * r_dim1;
	    z__2.r = c__[i__4] * r__[i__5].r, z__2.i = c__[i__4] * r__[i__5]
		    .i;
	    d_cnjg(&z__4, &s[ii]);
	    i__6 = i__ + j * r_dim1;
	    z__3.r = z__4.r * r__[i__6].r - z__4.i * r__[i__6].i, z__3.i = 
		    z__4.r * r__[i__6].i + z__4.i * r__[i__6].r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    r__[i__3].r = z__1.r, r__[i__3].i = z__1.i;
	    i__3 = i__ + j * r_dim1;
	    r__[i__3].r = t.r, r__[i__3].i = t.i;
/* L190: */
	}
L200:
	if (j >= *l) {
	    goto L210;
	}
	jj = j - *k + 1;
	i__2 = jj;
	t.r = s[i__2].r, t.i = s[i__2].i;
	crotg_(&r__[j + j * r_dim1], &t, &c__[jj], &s[jj]);
L210:
/* L220: */
	;
    }

/*        APPLY THE ROTATIONS TO Z. */

    if (*nz < 1) {
	goto L250;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = lm1;
	for (i__ = *k; i__ <= i__2; ++i__) {
	    ii = i__ - km1;
	    i__3 = ii;
	    i__4 = i__ + j * z_dim1;
	    z__2.r = c__[i__3] * z__[i__4].r, z__2.i = c__[i__3] * z__[i__4]
		    .i;
	    i__5 = ii;
	    i__6 = i__ + 1 + j * z_dim1;
	    z__3.r = s[i__5].r * z__[i__6].r - s[i__5].i * z__[i__6].i, 
		    z__3.i = s[i__5].r * z__[i__6].i + s[i__5].i * z__[i__6]
		    .r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = i__ + 1 + j * z_dim1;
	    i__4 = ii;
	    i__5 = i__ + 1 + j * z_dim1;
	    z__2.r = c__[i__4] * z__[i__5].r, z__2.i = c__[i__4] * z__[i__5]
		    .i;
	    d_cnjg(&z__4, &s[ii]);
	    i__6 = i__ + j * z_dim1;
	    z__3.r = z__4.r * z__[i__6].r - z__4.i * z__[i__6].i, z__3.i = 
		    z__4.r * z__[i__6].i + z__4.i * z__[i__6].r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = i__ + j * z_dim1;
	    z__[i__3].r = t.r, z__[i__3].i = t.i;
/* L230: */
	}
/* L240: */
    }
L250:
L260:
    return 0;
} /* cchex_ */

/* Subroutine */ int cqrdc_(doublecomplex *x, integer *ldx, integer *n, 
	integer *p, doublecomplex *qraux, integer *jpvt, doublecomplex *work, 
	integer *job)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *), z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), z_sqrt(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer j, l;
    static doublecomplex t;
    static integer jj, jp, pl, pu;
    static doublereal tt;
    static integer lp1, lup;
    static logical negj;
    static integer maxj;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static logical swapj;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublecomplex nrmxl;
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);
    static doublereal maxnrm;


/*     CQRDC USES HOUSEHOLDER TRANSFORMATIONS TO COMPUTE THE QR */
/*     FACTORIZATION OF AN N BY P MATRIX X.  COLUMN PIVOTING */
/*     BASED ON THE 2-NORMS OF THE REDUCED COLUMNS MAY BE */
/*     PERFORMED AT THE USERS OPTION. */

/*     ON ENTRY */

/*        X       COMPLEX(LDX,P), WHERE LDX .GE. N. */
/*                X CONTAINS THE MATRIX WHOSE DECOMPOSITION IS TO BE */
/*                COMPUTED. */

/*        LDX     INTEGER. */
/*                LDX IS THE LEADING DIMENSION OF THE ARRAY X. */

/*        N       INTEGER. */
/*                N IS THE NUMBER OF ROWS OF THE MATRIX X. */

/*        P       INTEGER. */
/*                P IS THE NUMBER OF COLUMNS OF THE MATRIX X. */

/*        JPVT    INTEGER(P). */
/*                JPVT CONTAINS INTEGERS THAT CONTROL THE SELECTION */
/*                OF THE PIVOT COLUMNS.  THE K-TH COLUMN X(K) OF X */
/*                IS PLACED IN ONE OF THREE CLASSES ACCORDING TO THE */
/*                VALUE OF JPVT(K). */

/*                   IF JPVT(K) .GT. 0, THEN X(K) IS AN INITIAL */
/*                                      COLUMN. */

/*                   IF JPVT(K) .EQ. 0, THEN X(K) IS A FREE COLUMN. */

/*                   IF JPVT(K) .LT. 0, THEN X(K) IS A FINAL COLUMN. */

/*                BEFORE THE DECOMPOSITION IS COMPUTED, INITIAL COLUMNS */
/*                ARE MOVED TO THE BEGINNING OF THE ARRAY X AND FINAL */
/*                COLUMNS TO THE END.  BOTH INITIAL AND FINAL COLUMNS */
/*                ARE FROZEN IN PLACE DURING THE COMPUTATION AND ONLY */
/*                FREE COLUMNS ARE MOVED.  AT THE K-TH STAGE OF THE */
/*                REDUCTION, IF X(K) IS OCCUPIED BY A FREE COLUMN */
/*                IT IS INTERCHANGED WITH THE FREE COLUMN OF LARGEST */
/*                REDUCED NORM.  JPVT IS NOT REFERENCED IF */
/*                JOB .EQ. 0. */

/*        WORK    COMPLEX(P). */
/*                WORK IS A WORK ARRAY.  WORK IS NOT REFERENCED IF */
/*                JOB .EQ. 0. */

/*        JOB     INTEGER. */
/*                JOB IS AN INTEGER THAT INITIATES COLUMN PIVOTING. */
/*                IF JOB .EQ. 0, NO PIVOTING IS DONE. */
/*                IF JOB .NE. 0, PIVOTING IS DONE. */

/*     ON RETURN */

/*        X       X CONTAINS IN ITS UPPER TRIANGLE THE UPPER */
/*                TRIANGULAR MATRIX R OF THE QR FACTORIZATION. */
/*                BELOW ITS DIAGONAL X CONTAINS INFORMATION FROM */
/*                WHICH THE UNITARY PART OF THE DECOMPOSITION */
/*                CAN BE RECOVERED.  NOTE THAT IF PIVOTING HAS */
/*                BEEN REQUESTED, THE DECOMPOSITION IS NOT THAT */
/*                OF THE ORIGINAL MATRIX X BUT THAT OF X */
/*                WITH ITS COLUMNS PERMUTED AS DESCRIBED BY JPVT. */

/*        QRAUX   COMPLEX(P). */
/*                QRAUX CONTAINS FURTHER INFORMATION REQUIRED TO RECOVER */
/*                THE UNITARY PART OF THE DECOMPOSITION. */

/*        JPVT    JPVT(K) CONTAINS THE INDEX OF THE COLUMN OF THE */
/*                ORIGINAL MATRIX THAT HAS BEEN INTERCHANGED INTO */
/*                THE K-TH COLUMN, IF PIVOTING WAS REQUESTED. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     CQRDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS. */

/*     BLAS CAXPY,CDOTC,CSCAL,CSWAP,SCNRM2 */
/*     FORTRAN ABS,DIMAG,MAX,CABS,DCMPLX,CSQRT,MIN0,DBLE */

/*     INTERNAL VARIABLES */



    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --qraux;
    --jpvt;
    --work;

    /* Function Body */
    pl = 1;
    pu = 0;
    if (*job == 0) {
	goto L60;
    }

/*        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS */
/*        ACCORDING TO JPVT. */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	swapj = jpvt[j] > 0;
	negj = jpvt[j] < 0;
	jpvt[j] = j;
	if (negj) {
	    jpvt[j] = -j;
	}
	if (! swapj) {
	    goto L10;
	}
	if (j != pl) {
	    cswap_(n, &x[pl * x_dim1 + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
	}
	jpvt[j] = jpvt[pl];
	jpvt[pl] = j;
	++pl;
L10:
/* L20: */
	;
    }
    pu = *p;
    i__1 = *p;
    for (jj = 1; jj <= i__1; ++jj) {
	j = *p - jj + 1;
	if (jpvt[j] >= 0) {
	    goto L40;
	}
	jpvt[j] = -jpvt[j];
	if (j == pu) {
	    goto L30;
	}
	cswap_(n, &x[pu * x_dim1 + 1], &c__1, &x[j * x_dim1 + 1], &c__1);
	jp = jpvt[pu];
	jpvt[pu] = jpvt[j];
	jpvt[j] = jp;
L30:
	--pu;
L40:
/* L50: */
	;
    }
L60:

/*     COMPUTE THE NORMS OF THE FREE COLUMNS. */

    if (pu < pl) {
	goto L80;
    }
    i__1 = pu;
    for (j = pl; j <= i__1; ++j) {
	i__2 = j;
	d__1 = scnrm2_(n, &x[j * x_dim1 + 1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	qraux[i__2].r = z__1.r, qraux[i__2].i = z__1.i;
	i__2 = j;
	i__3 = j;
	work[i__2].r = qraux[i__3].r, work[i__2].i = qraux[i__3].i;
/* L70: */
    }
L80:

/*     PERFORM THE HOUSEHOLDER REDUCTION OF X. */

    lup = min(*n,*p);
    i__1 = lup;
    for (l = 1; l <= i__1; ++l) {
	if (l < pl || l >= pu) {
	    goto L120;
	}

/*           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT */
/*           INTO THE PIVOT POSITION. */

	maxnrm = 0.;
	maxj = l;
	i__2 = pu;
	for (j = l; j <= i__2; ++j) {
	    i__3 = j;
	    if (qraux[i__3].r <= maxnrm) {
		goto L90;
	    }
	    i__3 = j;
	    maxnrm = qraux[i__3].r;
	    maxj = j;
L90:
/* L100: */
	    ;
	}
	if (maxj == l) {
	    goto L110;
	}
	cswap_(n, &x[l * x_dim1 + 1], &c__1, &x[maxj * x_dim1 + 1], &c__1);
	i__2 = maxj;
	i__3 = l;
	qraux[i__2].r = qraux[i__3].r, qraux[i__2].i = qraux[i__3].i;
	i__2 = maxj;
	i__3 = l;
	work[i__2].r = work[i__3].r, work[i__2].i = work[i__3].i;
	jp = jpvt[maxj];
	jpvt[maxj] = jpvt[l];
	jpvt[l] = jp;
L110:
L120:
	i__2 = l;
	qraux[i__2].r = 0., qraux[i__2].i = 0.;
	if (l == *n) {
	    goto L190;
	}

/*           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L. */

	i__2 = *n - l + 1;
	d__1 = scnrm2_(&i__2, &x[l + l * x_dim1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	nrmxl.r = z__1.r, nrmxl.i = z__1.i;
	if ((d__1 = nrmxl.r, abs(d__1)) + (d__2 = d_imag(&nrmxl), abs(d__2)) 
		== 0.) {
	    goto L180;
	}
	i__2 = l + l * x_dim1;
	if ((d__1 = x[i__2].r, abs(d__1)) + (d__2 = d_imag(&x[l + l * x_dim1])
		, abs(d__2)) != 0.) {
	    d__3 = z_abs(&nrmxl);
	    i__3 = l + l * x_dim1;
	    d__4 = z_abs(&x[l + l * x_dim1]);
	    z__2.r = x[i__3].r / d__4, z__2.i = x[i__3].i / d__4;
	    z__1.r = d__3 * z__2.r, z__1.i = d__3 * z__2.i;
	    nrmxl.r = z__1.r, nrmxl.i = z__1.i;
	}
	i__2 = *n - l + 1;
	z_div(&z__1, &c_b1092, &nrmxl);
	cscal_(&i__2, &z__1, &x[l + l * x_dim1], &c__1);
	i__2 = l + l * x_dim1;
	i__3 = l + l * x_dim1;
	z__1.r = x[i__3].r + 1., z__1.i = x[i__3].i + 0.;
	x[i__2].r = z__1.r, x[i__2].i = z__1.i;

/*              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS, */
/*              UPDATING THE NORMS. */

	lp1 = l + 1;
	if (*p < lp1) {
	    goto L170;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    cdotc_(&z__3, &i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1]
		    , &c__1);
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
	    z_div(&z__1, &z__2, &x[l + l * x_dim1]);
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = *n - l + 1;
	    caxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
	    if (j < pl || j > pu) {
		goto L150;
	    }
	    i__3 = j;
	    if ((d__1 = qraux[i__3].r, abs(d__1)) + (d__2 = d_imag(&qraux[j]),
		     abs(d__2)) == 0.) {
		goto L150;
	    }
	    i__3 = j;
/* Computing 2nd power */
	    d__1 = z_abs(&x[l + j * x_dim1]) / qraux[i__3].r;
	    tt = 1. - d__1 * d__1;
	    tt = max(tt,0.);
	    z__1.r = tt, z__1.i = 0.;
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = j;
	    i__4 = j;
/* Computing 2nd power */
	    d__1 = qraux[i__3].r / work[i__4].r;
	    tt = tt * .05f * (d__1 * d__1) + 1.;
	    if (tt == 1.) {
		goto L130;
	    }
	    i__3 = j;
	    i__4 = j;
	    z_sqrt(&z__2, &t);
	    z__1.r = qraux[i__4].r * z__2.r - qraux[i__4].i * z__2.i, z__1.i =
		     qraux[i__4].r * z__2.i + qraux[i__4].i * z__2.r;
	    qraux[i__3].r = z__1.r, qraux[i__3].i = z__1.i;
	    goto L140;
L130:
	    i__3 = j;
	    i__4 = *n - l;
	    d__1 = scnrm2_(&i__4, &x[l + 1 + j * x_dim1], &c__1);
	    z__1.r = d__1, z__1.i = 0.;
	    qraux[i__3].r = z__1.r, qraux[i__3].i = z__1.i;
	    i__3 = j;
	    i__4 = j;
	    work[i__3].r = qraux[i__4].r, work[i__3].i = qraux[i__4].i;
L140:
L150:
/* L160: */
	    ;
	}
L170:

/*              SAVE THE TRANSFORMATION. */

	i__2 = l;
	i__3 = l + l * x_dim1;
	qraux[i__2].r = x[i__3].r, qraux[i__2].i = x[i__3].i;
	i__2 = l + l * x_dim1;
	z__1.r = -nrmxl.r, z__1.i = -nrmxl.i;
	x[i__2].r = z__1.r, x[i__2].i = z__1.i;
L180:
L190:
/* L200: */
	;
    }
    return 0;
} /* cqrdc_ */

/* Subroutine */ int cqrsl_(doublecomplex *x, integer *ldx, integer *n, 
	integer *k, doublecomplex *qraux, doublecomplex *y, doublecomplex *qy,
	 doublecomplex *qty, doublecomplex *b, doublecomplex *rsd, 
	doublecomplex *xb, integer *job, integer *info)
{
    /* System generated locals */
    integer x_dim1, x_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, j;
    static doublecomplex t;
    static logical cb;
    static integer jj;
    static logical cr;
    static integer ju, kp1;
    static logical cxb, cqy;
    static doublecomplex temp;
    static logical cqty;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int ccopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     CQRSL APPLIES THE OUTPUT OF CQRDC TO COMPUTE COORDINATE */
/*     TRANSFORMATIONS, PROJECTIONS, AND LEAST SQUARES SOLUTIONS. */
/*     FOR K .LE. MIN(N,P), LET XK BE THE MATRIX */

/*            XK = (X(JPVT(1)),X(JPVT(2)), ... ,X(JPVT(K))) */

/*     FORMED FROM COLUMNNS JPVT(1), ... ,JPVT(K) OF THE ORIGINAL */
/*     N X P MATRIX X THAT WAS INPUT TO CQRDC (IF NO PIVOTING WAS */
/*     DONE, XK CONSISTS OF THE FIRST K COLUMNS OF X IN THEIR */
/*     ORIGINAL ORDER).  CQRDC PRODUCES A FACTORED UNITARY MATRIX Q */
/*     AND AN UPPER TRIANGULAR MATRIX R SUCH THAT */

/*              XK = Q * (R) */
/*                       (0) */

/*     THIS INFORMATION IS CONTAINED IN CODED FORM IN THE ARRAYS */
/*     X AND QRAUX. */

/*     ON ENTRY */

/*        X      COMPLEX(LDX,P). */
/*               X CONTAINS THE OUTPUT OF CQRDC. */

/*        LDX    INTEGER. */
/*               LDX IS THE LEADING DIMENSION OF THE ARRAY X. */

/*        N      INTEGER. */
/*               N IS THE NUMBER OF ROWS OF THE MATRIX XK.  IT MUST */
/*               HAVE THE SAME VALUE AS N IN CQRDC. */

/*        K      INTEGER. */
/*               K IS THE NUMBER OF COLUMNS OF THE MATRIX XK.  K */
/*               MUST NNOT BE GREATER THAN MIN(N,P), WHERE P IS THE */
/*               SAME AS IN THE CALLING SEQUENCE TO CQRDC. */

/*        QRAUX  COMPLEX(P). */
/*               QRAUX CONTAINS THE AUXILIARY OUTPUT FROM CQRDC. */

/*        Y      COMPLEX(N) */
/*               Y CONTAINS AN N-VECTOR THAT IS TO BE MANIPULATED */
/*               BY CQRSL. */

/*        JOB    INTEGER. */
/*               JOB SPECIFIES WHAT IS TO BE COMPUTED.  JOB HAS */
/*               THE DECIMAL EXPANSION ABCDE, WITH THE FOLLOWING */
/*               MEANING. */

/*                    IF A.NE.0, COMPUTE QY. */
/*                    IF B,C,D, OR E .NE. 0, COMPUTE QTY. */
/*                    IF C.NE.0, COMPUTE B. */
/*                    IF D.NE.0, COMPUTE RSD. */
/*                    IF E.NE.0, COMPUTE XB. */

/*               NOTE THAT A REQUEST TO COMPUTE B, RSD, OR XB */
/*               AUTOMATICALLY TRIGGERS THE COMPUTATION OF QTY, FOR */
/*               WHICH AN ARRAY MUST BE PROVIDED IN THE CALLING */
/*               SEQUENCE. */

/*     ON RETURN */

/*        QY     COMPLEX(N). */
/*               QY CONNTAINS Q*Y, IF ITS COMPUTATION HAS BEEN */
/*               REQUESTED. */

/*        QTY    COMPLEX(N). */
/*               QTY CONTAINS CTRANS(Q)*Y, IF ITS COMPUTATION HAS */
/*               BEEN REQUESTED.  HERE CTRANS(Q) IS THE CONJUGATE */
/*               TRANSPOSE OF THE MATRIX Q. */

/*        B      COMPLEX(K) */
/*               B CONTAINS THE SOLUTION OF THE LEAST SQUARES PROBLEM */

/*                    MINIMIZE NORM2(Y - XK*B), */

/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  (NOTE THAT */
/*               IF PIVOTING WAS REQUESTED IN CQRDC, THE J-TH */
/*               COMPONENT OF B WILL BE ASSOCIATED WITH COLUMN JPVT(J) */
/*               OF THE ORIGINAL MATRIX X THAT WAS INPUT INTO CQRDC.) */

/*        RSD    COMPLEX(N). */
/*               RSD CONTAINS THE LEAST SQUARES RESIDUAL Y - XK*B, */
/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  RSD IS */
/*               ALSO THE ORTHOGONAL PROJECTION OF Y ONTO THE */
/*               ORTHOGONAL COMPLEMENT OF THE COLUMN SPACE OF XK. */

/*        XB     COMPLEX(N). */
/*               XB CONTAINS THE LEAST SQUARES APPROXIMATION XK*B, */
/*               IF ITS COMPUTATION HAS BEEN REQUESTED.  XB IS ALSO */
/*               THE ORTHOGONAL PROJECTION OF Y ONTO THE COLUMN SPACE */
/*               OF X. */

/*        INFO   INTEGER. */
/*               INFO IS ZERO UNLESS THE COMPUTATION OF B HAS */
/*               BEEN REQUESTED AND R IS EXACTLY SINGULAR.  IN */
/*               THIS CASE, INFO IS THE INDEX OF THE FIRST ZERO */
/*               DIAGONAL ELEMENT OF R AND B IS LEFT UNALTERED. */

/*     THE PARAMETERS QY, QTY, B, RSD, AND XB ARE NOT REFERENCED */
/*     IF THEIR COMPUTATION IS NOT REQUESTED AND IN THIS CASE */
/*     CAN BE REPLACED BY DUMMY VARIABLES IN THE CALLING PROGRAM. */
/*     TO SAVE STORAGE, THE USER MAY IN SOME CASES USE THE SAME */
/*     ARRAY FOR DIFFERENT PARAMETERS IN THE CALLING SEQUENCE.  A */
/*     FREQUENTLY OCCURING EXAMPLE IS WHEN ONE WISHES TO COMPUTE */
/*     ANY OF B, RSD, OR XB AND DOES NOT NEED Y OR QTY.  IN THIS */
/*     CASE ONE MAY IDENTIFY Y, QTY, AND ONE OF B, RSD, OR XB, WHILE */
/*     PROVIDING SEPARATE ARRAYS FOR ANYTHING ELSE THAT IS TO BE */
/*     COMPUTED.  THUS THE CALLING SEQUENCE */

/*          CALL CQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO) */

/*     WILL RESULT IN THE COMPUTATION OF B AND RSD, WITH RSD */
/*     OVERWRITING Y.  MORE GENERALLY, EACH ITEM IN THE FOLLOWING */
/*     LIST CONTAINS GROUPS OF PERMISSIBLE IDENTIFICATIONS FOR */
/*     A SINGLE CALLINNG SEQUENCE. */

/*          1. (Y,QTY,B) (RSD) (XB) (QY) */

/*          2. (Y,QTY,RSD) (B) (XB) (QY) */

/*          3. (Y,QTY,XB) (B) (RSD) (QY) */

/*          4. (Y,QY) (QTY,B) (RSD) (XB) */

/*          5. (Y,QY) (QTY,RSD) (B) (XB) */

/*          6. (Y,QY) (QTY,XB) (B) (RSD) */

/*     IN ANY GROUP THE VALUE RETURNED IN THE ARRAY ALLOCATED TO */
/*     THE GROUP CORRESPONDS TO THE LAST MEMBER OF THE GROUP. */

/*     LINPACK. THIS VERSION DATED 08/14/78 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     CQRSL USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS. */

/*     BLAS CAXPY,CCOPY,CDOTC */
/*     FORTRAN ABS,DIMAG,MIN0,MOD,DBLE */

/*     INTERNAL VARIABLES */



/*     SET INFO FLAG. */

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --qraux;
    --y;
    --qy;
    --qty;
    --b;
    --rsd;
    --xb;

    /* Function Body */
    *info = 0;

/*     DETERMINE WHAT IS TO BE COMPUTED. */

    cqy = *job / 10000 != 0;
    cqty = *job % 10000 != 0;
    cb = *job % 1000 / 100 != 0;
    cr = *job % 100 / 10 != 0;
    cxb = *job % 10 != 0;
/* Computing MIN */
    i__1 = *k, i__2 = *n - 1;
    ju = min(i__1,i__2);

/*     SPECIAL ACTION WHEN N=1. */

    if (ju != 0) {
	goto L40;
    }
    if (cqy) {
	qy[1].r = y[1].r, qy[1].i = y[1].i;
    }
    if (cqty) {
	qty[1].r = y[1].r, qty[1].i = y[1].i;
    }
    if (cxb) {
	xb[1].r = y[1].r, xb[1].i = y[1].i;
    }
    if (! cb) {
	goto L30;
    }
    i__1 = x_dim1 + 1;
    if ((d__1 = x[i__1].r, abs(d__1)) + (d__2 = d_imag(&x[x_dim1 + 1]), abs(
	    d__2)) != 0.) {
	goto L10;
    }
    *info = 1;
    goto L20;
L10:
    z_div(&z__1, &y[1], &x[x_dim1 + 1]);
    b[1].r = z__1.r, b[1].i = z__1.i;
L20:
L30:
    if (cr) {
	rsd[1].r = 0., rsd[1].i = 0.;
    }
    goto L250;
L40:

/*        SET UP TO COMPUTE QY OR QTY. */

    if (cqy) {
	ccopy_(n, &y[1], &c__1, &qy[1], &c__1);
    }
    if (cqty) {
	ccopy_(n, &y[1], &c__1, &qty[1], &c__1);
    }
    if (! cqy) {
	goto L70;
    }

/*           COMPUTE QY. */

    i__1 = ju;
    for (jj = 1; jj <= i__1; ++jj) {
	j = ju - jj + 1;
	i__2 = j;
	if ((d__1 = qraux[i__2].r, abs(d__1)) + (d__2 = d_imag(&qraux[j]), 
		abs(d__2)) == 0.) {
	    goto L50;
	}
	i__2 = j + j * x_dim1;
	temp.r = x[i__2].r, temp.i = x[i__2].i;
	i__2 = j + j * x_dim1;
	i__3 = j;
	x[i__2].r = qraux[i__3].r, x[i__2].i = qraux[i__3].i;
	i__2 = *n - j + 1;
	cdotc_(&z__3, &i__2, &x[j + j * x_dim1], &c__1, &qy[j], &c__1);
	z__2.r = -z__3.r, z__2.i = -z__3.i;
	z_div(&z__1, &z__2, &x[j + j * x_dim1]);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = *n - j + 1;
	caxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qy[j], &c__1);
	i__2 = j + j * x_dim1;
	x[i__2].r = temp.r, x[i__2].i = temp.i;
L50:
/* L60: */
	;
    }
L70:
    if (! cqty) {
	goto L100;
    }

/*           COMPUTE CTRANS(Q)*Y. */

    i__1 = ju;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	if ((d__1 = qraux[i__2].r, abs(d__1)) + (d__2 = d_imag(&qraux[j]), 
		abs(d__2)) == 0.) {
	    goto L80;
	}
	i__2 = j + j * x_dim1;
	temp.r = x[i__2].r, temp.i = x[i__2].i;
	i__2 = j + j * x_dim1;
	i__3 = j;
	x[i__2].r = qraux[i__3].r, x[i__2].i = qraux[i__3].i;
	i__2 = *n - j + 1;
	cdotc_(&z__3, &i__2, &x[j + j * x_dim1], &c__1, &qty[j], &c__1);
	z__2.r = -z__3.r, z__2.i = -z__3.i;
	z_div(&z__1, &z__2, &x[j + j * x_dim1]);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = *n - j + 1;
	caxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &qty[j], &c__1);
	i__2 = j + j * x_dim1;
	x[i__2].r = temp.r, x[i__2].i = temp.i;
L80:
/* L90: */
	;
    }
L100:

/*        SET UP TO COMPUTE B, RSD, OR XB. */

    if (cb) {
	ccopy_(k, &qty[1], &c__1, &b[1], &c__1);
    }
    kp1 = *k + 1;
    if (cxb) {
	ccopy_(k, &qty[1], &c__1, &xb[1], &c__1);
    }
    if (cr && *k < *n) {
	i__1 = *n - *k;
	ccopy_(&i__1, &qty[kp1], &c__1, &rsd[kp1], &c__1);
    }
    if (! cxb || kp1 > *n) {
	goto L120;
    }
    i__1 = *n;
    for (i__ = kp1; i__ <= i__1; ++i__) {
	i__2 = i__;
	xb[i__2].r = 0., xb[i__2].i = 0.;
/* L110: */
    }
L120:
    if (! cr) {
	goto L140;
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	rsd[i__2].r = 0., rsd[i__2].i = 0.;
/* L130: */
    }
L140:
    if (! cb) {
	goto L190;
    }

/*           COMPUTE B. */

    i__1 = *k;
    for (jj = 1; jj <= i__1; ++jj) {
	j = *k - jj + 1;
	i__2 = j + j * x_dim1;
	if ((d__1 = x[i__2].r, abs(d__1)) + (d__2 = d_imag(&x[j + j * x_dim1])
		, abs(d__2)) != 0.) {
	    goto L150;
	}
	*info = j;
/*           ......EXIT */
	goto L180;
L150:
	i__2 = j;
	z_div(&z__1, &b[j], &x[j + j * x_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	if (j == 1) {
	    goto L160;
	}
	i__2 = j;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &t, &x[j * x_dim1 + 1], &c__1, &b[1], &c__1);
L160:
/* L170: */
	;
    }
L180:
L190:
    if (! cr && ! cxb) {
	goto L240;
    }

/*           COMPUTE RSD OR XB AS REQUIRED. */

    i__1 = ju;
    for (jj = 1; jj <= i__1; ++jj) {
	j = ju - jj + 1;
	i__2 = j;
	if ((d__1 = qraux[i__2].r, abs(d__1)) + (d__2 = d_imag(&qraux[j]), 
		abs(d__2)) == 0.) {
	    goto L220;
	}
	i__2 = j + j * x_dim1;
	temp.r = x[i__2].r, temp.i = x[i__2].i;
	i__2 = j + j * x_dim1;
	i__3 = j;
	x[i__2].r = qraux[i__3].r, x[i__2].i = qraux[i__3].i;
	if (! cr) {
	    goto L200;
	}
	i__2 = *n - j + 1;
	cdotc_(&z__3, &i__2, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1);
	z__2.r = -z__3.r, z__2.i = -z__3.i;
	z_div(&z__1, &z__2, &x[j + j * x_dim1]);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = *n - j + 1;
	caxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &rsd[j], &c__1);
L200:
	if (! cxb) {
	    goto L210;
	}
	i__2 = *n - j + 1;
	cdotc_(&z__3, &i__2, &x[j + j * x_dim1], &c__1, &xb[j], &c__1);
	z__2.r = -z__3.r, z__2.i = -z__3.i;
	z_div(&z__1, &z__2, &x[j + j * x_dim1]);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = *n - j + 1;
	caxpy_(&i__2, &t, &x[j + j * x_dim1], &c__1, &xb[j], &c__1);
L210:
	i__2 = j + j * x_dim1;
	x[i__2].r = temp.r, x[i__2].i = temp.i;
L220:
/* L230: */
	;
    }
L240:
L250:
    return 0;
} /* cqrsl_ */

/* Subroutine */ int csvdc_(doublecomplex *x, integer *ldx, integer *n, 
	integer *p, doublecomplex *s, doublecomplex *e, doublecomplex *u, 
	integer *ldu, doublecomplex *v, integer *ldv, doublecomplex *work, 
	integer *job, integer *info)
{
    /* System generated locals */
    integer x_dim1, x_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double d_imag(doublecomplex *), z_abs(doublecomplex *);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);
    double sqrt(doublereal);

    /* Local variables */
    static doublereal b, c__, f, g;
    static integer i__, j, k, l, m;
    static doublecomplex r__, t;
    static doublereal t1, el;
    static integer kk;
    static doublereal cs;
    static integer ll, mm, ls;
    static doublereal sl;
    static integer lu;
    static doublereal sm, sn;
    static integer lm1, mm1, lp1, mp1, nct, ncu, lls, nrt;
    static doublereal emm1, smm1;
    static integer kase, jobu, iter;
    static doublereal test;
    static integer nctp1, nrtp1;
    extern /* Subroutine */ int cscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *);
    static doublereal scale;
    extern /* Double Complex */ VOID cdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    static doublereal shift;
    extern /* Subroutine */ int cswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    static integer maxit;
    extern /* Subroutine */ int caxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), csrot_(
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, doublereal *);
    static logical wantu, wantv;
    extern /* Subroutine */ int srotg_(doublereal *, doublereal *, doublereal 
	    *, doublereal *);
    static doublereal ztest;
    extern doublereal scnrm2_(integer *, doublecomplex *, integer *);



/*     CSVDC IS A SUBROUTINE TO REDUCE A COMPLEX NXP MATRIX X BY */
/*     UNITARY TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE */
/*     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE */
/*     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS, */
/*     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS. */

/*     ON ENTRY */

/*         X         COMPLEX(LDX,P), WHERE LDX.GE.N. */
/*                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE */
/*                   DECOMPOSITION IS TO BE COMPUTED.  X IS */
/*                   DESTROYED BY CSVDC. */

/*         LDX       INTEGER. */
/*                   LDX IS THE LEADING DIMENSION OF THE ARRAY X. */

/*         N         INTEGER. */
/*                   N IS THE NUMBER OF COLUMNS OF THE MATRIX X. */

/*         P         INTEGER. */
/*                   P IS THE NUMBER OF ROWS OF THE MATRIX X. */

/*         LDU       INTEGER. */
/*                   LDU IS THE LEADING DIMENSION OF THE ARRAY U */
/*                   (SEE BELOW). */

/*         LDV       INTEGER. */
/*                   LDV IS THE LEADING DIMENSION OF THE ARRAY V */
/*                   (SEE BELOW). */

/*         WORK      COMPLEX(N). */
/*                   WORK IS A SCRATCH ARRAY. */

/*         JOB       INTEGER. */
/*                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR */
/*                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB */
/*                   WITH THE FOLLOWING MEANING */

/*                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR */
/*                                  VECTORS. */
/*                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS */
/*                                  IN U. */
/*                        A.GE.2    RETURNS THE FIRST MIN(N,P) */
/*                                  LEFT SINGULAR VECTORS IN U. */
/*                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR */
/*                                  VECTORS. */
/*                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS */
/*                                  IN V. */

/*     ON RETURN */

/*         S         COMPLEX(MM), WHERE MM=MIN(N+1,P). */
/*                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE */
/*                   SINGULAR VALUES OF X ARRANGED IN DESCENDING */
/*                   ORDER OF MAGNITUDE. */

/*         E         COMPLEX(P). */
/*                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE */
/*                   DISCUSSION OF INFO FOR EXCEPTIONS. */

/*         U         COMPLEX(LDU,K), WHERE LDU.GE.N.  IF JOBA.EQ.1 THEN */
/*                                   K.EQ.N, IF JOBA.GE.2 THEN */
/*                                   K.EQ.MIN(N,P). */
/*                   U CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS. */
/*                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P */
/*                   OR IF JOBA.GT.2, THEN U MAY BE IDENTIFIED WITH X */
/*                   IN THE SUBROUTINE CALL. */

/*         V         COMPLEX(LDV,P), WHERE LDV.GE.P. */
/*                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS. */
/*                   V IS NOT REFERENCED IF JOBB.EQ.0.  IF P.LE.N, */
/*                   THEN V MAY BE IDENTIFIED WHTH X IN THE */
/*                   SUBROUTINE CALL. */

/*         INFO      INTEGER. */
/*                   THE SINGULAR VALUES (AND THEIR CORRESPONDING */
/*                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M) */
/*                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF */
/*                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR */
/*                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX */
/*                   B = CTRANS(U)*X*V IS THE BIDIAGONAL MATRIX */
/*                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE */
/*                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (CTRANS(U) */
/*                   IS THE CONJUGATE-TRANSPOSE OF U).  THUS THE */
/*                   SINGULAR VALUES OF X AND B ARE THE SAME. */

/*     LINPACK. THIS VERSION DATED 03/19/79 . */
/*     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB. */

/*     CSVDC USES THE FOLLOWING FUNCTIONS AND SUBPROGRAMS. */

/*     EXTERNAL CSROT */
/*     BLAS CAXPY,CDOTC,CSCAL,CSWAP,SCNRM2,SROTG */
/*     FORTRAN ABS,DIMAG,MAX,CABS,DCMPLX */
/*     FORTRAN CONJG,MAX0,MIN0,MOD,DBLE,SQRT */

/*     INTERNAL VARIABLES */



/*     SET THE MAXIMUM NUMBER OF ITERATIONS. */

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --s;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --work;

    /* Function Body */
    maxit = 30;

/*     DETERMINE WHAT IS TO BE COMPUTED. */

    wantu = FALSE_;
    wantv = FALSE_;
    jobu = *job % 100 / 10;
    ncu = *n;
    if (jobu > 1) {
	ncu = min(*n,*p);
    }
    if (jobu != 0) {
	wantu = TRUE_;
    }
    if (*job % 10 != 0) {
	wantv = TRUE_;
    }

/*     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS */
/*     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E. */

    *info = 0;
/* Computing MIN */
    i__1 = *n - 1;
    nct = min(i__1,*p);
/* Computing MAX */
/* Computing MIN */
    i__3 = *p - 2;
    i__1 = 0, i__2 = min(i__3,*n);
    nrt = max(i__1,i__2);
    lu = max(nct,nrt);
    if (lu < 1) {
	goto L170;
    }
    i__1 = lu;
    for (l = 1; l <= i__1; ++l) {
	lp1 = l + 1;
	if (l > nct) {
	    goto L20;
	}

/*           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND */
/*           PLACE THE L-TH DIAGONAL IN S(L). */

	i__2 = l;
	i__3 = *n - l + 1;
	d__1 = scnrm2_(&i__3, &x[l + l * x_dim1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
	i__2 = l;
	if ((d__1 = s[i__2].r, abs(d__1)) + (d__2 = d_imag(&s[l]), abs(d__2)) 
		== 0.) {
	    goto L10;
	}
	i__2 = l + l * x_dim1;
	if ((d__1 = x[i__2].r, abs(d__1)) + (d__2 = d_imag(&x[l + l * x_dim1])
		, abs(d__2)) != 0.) {
	    i__3 = l;
	    d__3 = z_abs(&s[l]);
	    i__4 = l + l * x_dim1;
	    d__4 = z_abs(&x[l + l * x_dim1]);
	    z__2.r = x[i__4].r / d__4, z__2.i = x[i__4].i / d__4;
	    z__1.r = d__3 * z__2.r, z__1.i = d__3 * z__2.i;
	    s[i__3].r = z__1.r, s[i__3].i = z__1.i;
	}
	i__2 = *n - l + 1;
	z_div(&z__1, &c_b1092, &s[l]);
	cscal_(&i__2, &z__1, &x[l + l * x_dim1], &c__1);
	i__2 = l + l * x_dim1;
	i__3 = l + l * x_dim1;
	z__1.r = x[i__3].r + 1., z__1.i = x[i__3].i + 0.;
	x[i__2].r = z__1.r, x[i__2].i = z__1.i;
L10:
	i__2 = l;
	i__3 = l;
	z__1.r = -s[i__3].r, z__1.i = -s[i__3].i;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
L20:
	if (*p < lp1) {
	    goto L50;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    if (l > nct) {
		goto L30;
	    }
	    i__3 = l;
	    if ((d__1 = s[i__3].r, abs(d__1)) + (d__2 = d_imag(&s[l]), abs(
		    d__2)) == 0.) {
		goto L30;
	    }

/*              APPLY THE TRANSFORMATION. */

	    i__3 = *n - l + 1;
	    cdotc_(&z__3, &i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1]
		    , &c__1);
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
	    z_div(&z__1, &z__2, &x[l + l * x_dim1]);
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = *n - l + 1;
	    caxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
L30:

/*           PLACE THE L-TH ROW OF X INTO  E FOR THE */
/*           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION. */

	    i__3 = j;
	    d_cnjg(&z__1, &x[l + j * x_dim1]);
	    e[i__3].r = z__1.r, e[i__3].i = z__1.i;
/* L40: */
	}
L50:
	if (! wantu || l > nct) {
	    goto L70;
	}

/*           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK */
/*           MULTIPLICATION. */

	i__2 = *n;
	for (i__ = l; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * u_dim1;
	    i__4 = i__ + l * x_dim1;
	    u[i__3].r = x[i__4].r, u[i__3].i = x[i__4].i;
/* L60: */
	}
L70:
	if (l > nrt) {
	    goto L150;
	}

/*           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE */
/*           L-TH SUPER-DIAGONAL IN E(L). */

	i__2 = l;
	i__3 = *p - l;
	d__1 = scnrm2_(&i__3, &e[lp1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
	i__2 = l;
	if ((d__1 = e[i__2].r, abs(d__1)) + (d__2 = d_imag(&e[l]), abs(d__2)) 
		== 0.) {
	    goto L80;
	}
	i__2 = lp1;
	if ((d__1 = e[i__2].r, abs(d__1)) + (d__2 = d_imag(&e[lp1]), abs(d__2)
		) != 0.) {
	    i__3 = l;
	    d__3 = z_abs(&e[l]);
	    i__4 = lp1;
	    d__4 = z_abs(&e[lp1]);
	    z__2.r = e[i__4].r / d__4, z__2.i = e[i__4].i / d__4;
	    z__1.r = d__3 * z__2.r, z__1.i = d__3 * z__2.i;
	    e[i__3].r = z__1.r, e[i__3].i = z__1.i;
	}
	i__2 = *p - l;
	z_div(&z__1, &c_b1092, &e[l]);
	cscal_(&i__2, &z__1, &e[lp1], &c__1);
	i__2 = lp1;
	i__3 = lp1;
	z__1.r = e[i__3].r + 1., z__1.i = e[i__3].i + 0.;
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
L80:
	i__2 = l;
	d_cnjg(&z__2, &e[l]);
	z__1.r = -z__2.r, z__1.i = -z__2.i;
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
	i__2 = l;
	if (lp1 > *n || (d__1 = e[i__2].r, abs(d__1)) + (d__2 = d_imag(&e[l]),
		 abs(d__2)) == 0.) {
	    goto L120;
	}

/*              APPLY THE TRANSFORMATION. */

	i__2 = *n;
	for (i__ = lp1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    work[i__3].r = 0., work[i__3].i = 0.;
/* L90: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    caxpy_(&i__3, &e[j], &x[lp1 + j * x_dim1], &c__1, &work[lp1], &
		    c__1);
/* L100: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    i__4 = j;
	    z__3.r = -e[i__4].r, z__3.i = -e[i__4].i;
	    z_div(&z__2, &z__3, &e[lp1]);
	    d_cnjg(&z__1, &z__2);
	    caxpy_(&i__3, &z__1, &work[lp1], &c__1, &x[lp1 + j * x_dim1], &
		    c__1);
/* L110: */
	}
L120:
	if (! wantv) {
	    goto L140;
	}

/*              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT */
/*              BACK MULTIPLICATION. */

	i__2 = *p;
	for (i__ = lp1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * v_dim1;
	    i__4 = i__;
	    v[i__3].r = e[i__4].r, v[i__3].i = e[i__4].i;
/* L130: */
	}
L140:
L150:
/* L160: */
	;
    }
L170:

/*     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M. */

/* Computing MIN */
    i__1 = *p, i__2 = *n + 1;
    m = min(i__1,i__2);
    nctp1 = nct + 1;
    nrtp1 = nrt + 1;
    if (nct < *p) {
	i__1 = nctp1;
	i__2 = nctp1 + nctp1 * x_dim1;
	s[i__1].r = x[i__2].r, s[i__1].i = x[i__2].i;
    }
    if (*n < m) {
	i__1 = m;
	s[i__1].r = 0., s[i__1].i = 0.;
    }
    if (nrtp1 < m) {
	i__1 = nrtp1;
	i__2 = nrtp1 + m * x_dim1;
	e[i__1].r = x[i__2].r, e[i__1].i = x[i__2].i;
    }
    i__1 = m;
    e[i__1].r = 0., e[i__1].i = 0.;

/*     IF REQUIRED, GENERATE U. */

    if (! wantu) {
	goto L300;
    }
    if (ncu < nctp1) {
	goto L200;
    }
    i__1 = ncu;
    for (j = nctp1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * u_dim1;
	    u[i__3].r = 0., u[i__3].i = 0.;
/* L180: */
	}
	i__2 = j + j * u_dim1;
	u[i__2].r = 1., u[i__2].i = 0.;
/* L190: */
    }
L200:
    if (nct < 1) {
	goto L290;
    }
    i__1 = nct;
    for (ll = 1; ll <= i__1; ++ll) {
	l = nct - ll + 1;
	i__2 = l;
	if ((d__1 = s[i__2].r, abs(d__1)) + (d__2 = d_imag(&s[l]), abs(d__2)) 
		== 0.) {
	    goto L250;
	}
	lp1 = l + 1;
	if (ncu < lp1) {
	    goto L220;
	}
	i__2 = ncu;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    cdotc_(&z__3, &i__3, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1]
		    , &c__1);
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
	    z_div(&z__1, &z__2, &u[l + l * u_dim1]);
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = *n - l + 1;
	    caxpy_(&i__3, &t, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &
		    c__1);
/* L210: */
	}
L220:
	i__2 = *n - l + 1;
	cscal_(&i__2, &c_b2561, &u[l + l * u_dim1], &c__1);
	i__2 = l + l * u_dim1;
	i__3 = l + l * u_dim1;
	z__1.r = u[i__3].r + 1., z__1.i = u[i__3].i + 0.;
	u[i__2].r = z__1.r, u[i__2].i = z__1.i;
	lm1 = l - 1;
	if (lm1 < 1) {
	    goto L240;
	}
	i__2 = lm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * u_dim1;
	    u[i__3].r = 0., u[i__3].i = 0.;
/* L230: */
	}
L240:
	goto L270;
L250:
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * u_dim1;
	    u[i__3].r = 0., u[i__3].i = 0.;
/* L260: */
	}
	i__2 = l + l * u_dim1;
	u[i__2].r = 1., u[i__2].i = 0.;
L270:
/* L280: */
	;
    }
L290:
L300:

/*     IF IT IS REQUIRED, GENERATE V. */

    if (! wantv) {
	goto L350;
    }
    i__1 = *p;
    for (ll = 1; ll <= i__1; ++ll) {
	l = *p - ll + 1;
	lp1 = l + 1;
	if (l > nrt) {
	    goto L320;
	}
	i__2 = l;
	if ((d__1 = e[i__2].r, abs(d__1)) + (d__2 = d_imag(&e[l]), abs(d__2)) 
		== 0.) {
	    goto L320;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *p - l;
	    cdotc_(&z__3, &i__3, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1);
	    z__2.r = -z__3.r, z__2.i = -z__3.i;
	    z_div(&z__1, &z__2, &v[lp1 + l * v_dim1]);
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = *p - l;
	    caxpy_(&i__3, &t, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1);
/* L310: */
	}
L320:
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * v_dim1;
	    v[i__3].r = 0., v[i__3].i = 0.;
/* L330: */
	}
	i__2 = l + l * v_dim1;
	v[i__2].r = 1., v[i__2].i = 0.;
/* L340: */
    }
L350:

/*     TRANSFORM S AND E SO THAT THEY ARE REAL. */

    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	if ((d__1 = s[i__2].r, abs(d__1)) + (d__2 = d_imag(&s[i__]), abs(d__2)
		) == 0.) {
	    goto L360;
	}
	d__1 = z_abs(&s[i__]);
	z__1.r = d__1, z__1.i = 0.;
	t.r = z__1.r, t.i = z__1.i;
	z_div(&z__1, &s[i__], &t);
	r__.r = z__1.r, r__.i = z__1.i;
	i__2 = i__;
	s[i__2].r = t.r, s[i__2].i = t.i;
	if (i__ < m) {
	    i__2 = i__;
	    z_div(&z__1, &e[i__], &r__);
	    e[i__2].r = z__1.r, e[i__2].i = z__1.i;
	}
	if (wantu) {
	    cscal_(n, &r__, &u[i__ * u_dim1 + 1], &c__1);
	}
L360:
/*     ...EXIT */
	if (i__ == m) {
	    goto L390;
	}
	i__2 = i__;
	if ((d__1 = e[i__2].r, abs(d__1)) + (d__2 = d_imag(&e[i__]), abs(d__2)
		) == 0.) {
	    goto L370;
	}
	d__1 = z_abs(&e[i__]);
	z__1.r = d__1, z__1.i = 0.;
	t.r = z__1.r, t.i = z__1.i;
	z_div(&z__1, &t, &e[i__]);
	r__.r = z__1.r, r__.i = z__1.i;
	i__2 = i__;
	e[i__2].r = t.r, e[i__2].i = t.i;
	i__2 = i__ + 1;
	i__3 = i__ + 1;
	z__1.r = s[i__3].r * r__.r - s[i__3].i * r__.i, z__1.i = s[i__3].r * 
		r__.i + s[i__3].i * r__.r;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
	if (wantv) {
	    cscal_(p, &r__, &v[(i__ + 1) * v_dim1 + 1], &c__1);
	}
L370:
/* L380: */
	;
    }
L390:

/*     MAIN ITERATION LOOP FOR THE SINGULAR VALUES. */

    mm = m;
    iter = 0;
L400:

/*        QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND. */

/*     ...EXIT */
    if (m == 0) {
	goto L660;
    }

/*        IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET */
/*        FLAG AND RETURN. */

    if (iter < maxit) {
	goto L410;
    }
    *info = m;
/*     ......EXIT */
    goto L660;
L410:

/*        THIS SECTION OF THE PROGRAM INSPECTS FOR */
/*        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON */
/*        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS. */

/*           KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M */
/*           KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M */
/*           KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND */
/*                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP). */
/*           KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE). */

    i__1 = m;
    for (ll = 1; ll <= i__1; ++ll) {
	l = m - ll;
/*        ...EXIT */
	if (l == 0) {
	    goto L440;
	}
	test = z_abs(&s[l]) + z_abs(&s[l + 1]);
	ztest = test + z_abs(&e[l]);
	if (ztest != test) {
	    goto L420;
	}
	i__2 = l;
	e[i__2].r = 0., e[i__2].i = 0.;
/*        ......EXIT */
	goto L440;
L420:
/* L430: */
	;
    }
L440:
    if (l != m - 1) {
	goto L450;
    }
    kase = 4;
    goto L520;
L450:
    lp1 = l + 1;
    mp1 = m + 1;
    i__1 = mp1;
    for (lls = lp1; lls <= i__1; ++lls) {
	ls = m - lls + lp1;
/*           ...EXIT */
	if (ls == l) {
	    goto L480;
	}
	test = 0.;
	if (ls != m) {
	    test += z_abs(&e[ls]);
	}
	if (ls != l + 1) {
	    test += z_abs(&e[ls - 1]);
	}
	ztest = test + z_abs(&s[ls]);
	if (ztest != test) {
	    goto L460;
	}
	i__2 = ls;
	s[i__2].r = 0., s[i__2].i = 0.;
/*           ......EXIT */
	goto L480;
L460:
/* L470: */
	;
    }
L480:
    if (ls != l) {
	goto L490;
    }
    kase = 3;
    goto L510;
L490:
    if (ls != m) {
	goto L500;
    }
    kase = 1;
    goto L510;
L500:
    kase = 2;
    l = ls;
L510:
L520:
    ++l;

/*        PERFORM THE TASK INDICATED BY KASE. */

    switch (kase) {
	case 1:  goto L530;
	case 2:  goto L560;
	case 3:  goto L580;
	case 4:  goto L610;
    }

/*        DEFLATE NEGLIGIBLE S(M). */

L530:
    mm1 = m - 1;
    i__1 = m - 1;
    f = e[i__1].r;
    i__1 = m - 1;
    e[i__1].r = 0., e[i__1].i = 0.;
    i__1 = mm1;
    for (kk = l; kk <= i__1; ++kk) {
	k = mm1 - kk + l;
	i__2 = k;
	t1 = s[i__2].r;
	srotg_(&t1, &f, &cs, &sn);
	i__2 = k;
	z__1.r = t1, z__1.i = 0.;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
	if (k == l) {
	    goto L540;
	}
	i__2 = k - 1;
	f = -sn * e[i__2].r;
	i__2 = k - 1;
	i__3 = k - 1;
	z__1.r = cs * e[i__3].r, z__1.i = cs * e[i__3].i;
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
L540:
	if (wantv) {
	    csrot_(p, &v[k * v_dim1 + 1], &c__1, &v[m * v_dim1 + 1], &c__1, &
		    cs, &sn);
	}
/* L550: */
    }
    goto L650;

/*        SPLIT AT NEGLIGIBLE S(L). */

L560:
    i__1 = l - 1;
    f = e[i__1].r;
    i__1 = l - 1;
    e[i__1].r = 0., e[i__1].i = 0.;
    i__1 = m;
    for (k = l; k <= i__1; ++k) {
	i__2 = k;
	t1 = s[i__2].r;
	srotg_(&t1, &f, &cs, &sn);
	i__2 = k;
	z__1.r = t1, z__1.i = 0.;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
	i__2 = k;
	f = -sn * e[i__2].r;
	i__2 = k;
	i__3 = k;
	z__1.r = cs * e[i__3].r, z__1.i = cs * e[i__3].i;
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
	if (wantu) {
	    csrot_(n, &u[k * u_dim1 + 1], &c__1, &u[(l - 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L570: */
    }
    goto L650;

/*        PERFORM ONE QR STEP. */

L580:

/*           CALCULATE THE SHIFT. */

/* Computing MAX */
    d__1 = z_abs(&s[m]), d__2 = z_abs(&s[m - 1]), d__1 = max(d__1,d__2), d__2 
	    = z_abs(&e[m - 1]), d__1 = max(d__1,d__2), d__2 = z_abs(&s[l]), 
	    d__1 = max(d__1,d__2), d__2 = z_abs(&e[l]);
    scale = max(d__1,d__2);
    i__1 = m;
    sm = s[i__1].r / scale;
    i__1 = m - 1;
    smm1 = s[i__1].r / scale;
    i__1 = m - 1;
    emm1 = e[i__1].r / scale;
    i__1 = l;
    sl = s[i__1].r / scale;
    i__1 = l;
    el = e[i__1].r / scale;
/* Computing 2nd power */
    d__1 = emm1;
    b = ((smm1 + sm) * (smm1 - sm) + d__1 * d__1) / 2.;
/* Computing 2nd power */
    d__1 = sm * emm1;
    c__ = d__1 * d__1;
    shift = 0.;
    if (b == 0. && c__ == 0.) {
	goto L590;
    }
/* Computing 2nd power */
    d__1 = b;
    shift = sqrt(d__1 * d__1 + c__);
    if (b < 0.) {
	shift = -shift;
    }
    shift = c__ / (b + shift);
L590:
    f = (sl + sm) * (sl - sm) - shift;
    g = sl * el;

/*           CHASE ZEROS. */

    mm1 = m - 1;
    i__1 = mm1;
    for (k = l; k <= i__1; ++k) {
	srotg_(&f, &g, &cs, &sn);
	if (k != l) {
	    i__2 = k - 1;
	    z__1.r = f, z__1.i = 0.;
	    e[i__2].r = z__1.r, e[i__2].i = z__1.i;
	}
	i__2 = k;
	i__3 = k;
	f = cs * s[i__2].r + sn * e[i__3].r;
	i__2 = k;
	i__3 = k;
	z__2.r = cs * e[i__3].r, z__2.i = cs * e[i__3].i;
	i__4 = k;
	z__3.r = sn * s[i__4].r, z__3.i = sn * s[i__4].i;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
	i__2 = k + 1;
	g = sn * s[i__2].r;
	i__2 = k + 1;
	i__3 = k + 1;
	z__1.r = cs * s[i__3].r, z__1.i = cs * s[i__3].i;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
	if (wantv) {
	    csrot_(p, &v[k * v_dim1 + 1], &c__1, &v[(k + 1) * v_dim1 + 1], &
		    c__1, &cs, &sn);
	}
	srotg_(&f, &g, &cs, &sn);
	i__2 = k;
	z__1.r = f, z__1.i = 0.;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
	i__2 = k;
	i__3 = k + 1;
	f = cs * e[i__2].r + sn * s[i__3].r;
	i__2 = k + 1;
	d__1 = -sn;
	i__3 = k;
	z__2.r = d__1 * e[i__3].r, z__2.i = d__1 * e[i__3].i;
	i__4 = k + 1;
	z__3.r = cs * s[i__4].r, z__3.i = cs * s[i__4].i;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
	i__2 = k + 1;
	g = sn * e[i__2].r;
	i__2 = k + 1;
	i__3 = k + 1;
	z__1.r = cs * e[i__3].r, z__1.i = cs * e[i__3].i;
	e[i__2].r = z__1.r, e[i__2].i = z__1.i;
	if (wantu && k < *n) {
	    csrot_(n, &u[k * u_dim1 + 1], &c__1, &u[(k + 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L600: */
    }
    i__1 = m - 1;
    z__1.r = f, z__1.i = 0.;
    e[i__1].r = z__1.r, e[i__1].i = z__1.i;
    ++iter;
    goto L650;

/*        CONVERGENCE. */

L610:

/*           MAKE THE SINGULAR VALUE  POSITIVE */

    i__1 = l;
    if (s[i__1].r >= 0.) {
	goto L620;
    }
    i__1 = l;
    i__2 = l;
    z__1.r = -s[i__2].r, z__1.i = -s[i__2].i;
    s[i__1].r = z__1.r, s[i__1].i = z__1.i;
    if (wantv) {
	cscal_(p, &c_b2561, &v[l * v_dim1 + 1], &c__1);
    }
L620:

/*           ORDER THE SINGULAR VALUE. */

L630:
    if (l == mm) {
	goto L640;
    }
/*           ...EXIT */
    i__1 = l;
    i__2 = l + 1;
    if (s[i__1].r >= s[i__2].r) {
	goto L640;
    }
    i__1 = l;
    t.r = s[i__1].r, t.i = s[i__1].i;
    i__1 = l;
    i__2 = l + 1;
    s[i__1].r = s[i__2].r, s[i__1].i = s[i__2].i;
    i__1 = l + 1;
    s[i__1].r = t.r, s[i__1].i = t.i;
    if (wantv && l < *p) {
	cswap_(p, &v[l * v_dim1 + 1], &c__1, &v[(l + 1) * v_dim1 + 1], &c__1);
    }
    if (wantu && l < *n) {
	cswap_(n, &u[l * u_dim1 + 1], &c__1, &u[(l + 1) * u_dim1 + 1], &c__1);
    }
    ++l;
    goto L630;
L640:
    iter = 0;
    --m;
L650:
    goto L400;
L660:
    return 0;
} /* csvdc_ */

/* Subroutine */ int csrot_(integer *n, doublecomplex *cx, integer *incx, 
	doublecomplex *cy, integer *incy, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublecomplex z__1, z__2, z__3;

    /* Local variables */
    static integer i__, ix, iy;
    static doublecomplex ctemp;


/*     APPLIES A PLANE ROTATION, WHERE THE COS AND SIN (C AND S) ARE REAL */
/*     AND THE VECTORS CX AND CY ARE COMPLEX. */
/*     JACK DONGARRA, LINPACK, 3/11/78. */


    /* Parameter adjustments */
    --cy;
    --cx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL */
/*         TO 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = ix;
	z__2.r = *c__ * cx[i__2].r, z__2.i = *c__ * cx[i__2].i;
	i__3 = iy;
	z__3.r = *s * cy[i__3].r, z__3.i = *s * cy[i__3].i;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	i__2 = iy;
	i__3 = iy;
	z__2.r = *c__ * cy[i__3].r, z__2.i = *c__ * cy[i__3].i;
	i__4 = ix;
	z__3.r = *s * cx[i__4].r, z__3.i = *s * cx[i__4].i;
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
	z__2.r = *c__ * cx[i__2].r, z__2.i = *c__ * cx[i__2].i;
	i__3 = i__;
	z__3.r = *s * cy[i__3].r, z__3.i = *s * cy[i__3].i;
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	ctemp.r = z__1.r, ctemp.i = z__1.i;
	i__2 = i__;
	i__3 = i__;
	z__2.r = *c__ * cy[i__3].r, z__2.i = *c__ * cy[i__3].i;
	i__4 = i__;
	z__3.r = *s * cx[i__4].r, z__3.i = *s * cx[i__4].i;
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	cy[i__2].r = z__1.r, cy[i__2].i = z__1.i;
	i__2 = i__;
	cx[i__2].r = ctemp.r, cx[i__2].i = ctemp.i;
/* L30: */
    }
    return 0;
} /* csrot_ */

