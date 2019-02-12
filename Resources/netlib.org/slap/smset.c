/* smset.f -- translated by f2c (version 20100827).
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

/* Common Block Declarations */

struct {
    real soln[1];
} solblk_;

#define solblk_1 solblk_

/* Table of constant values */

static integer c__53 = 53;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__0 = 0;
static real c_b44 = 0.f;

/* DECK SSDS */
/* Subroutine */ int ssds_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, real *dinv)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer icol;

/* ***BEGIN PROLOGUE  SSDS */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSDS-S), */
/*             SLAP Sparse, Diagonal */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Diagonal Scaling Preconditioner SLAP Set Up. */
/*            Routine to compute the inverse of the diagonal of a matrix */
/*            stored in the SLAP Column format. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM */
/*     REAL    A(NELT), DINV(N) */

/*     CALL SSDS( N, NELT, IA, JA, A, ISYM, DINV ) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of elements in arrays IA, JA, and A. */
/* IA     :INOUT    Integer IA(NELT). */
/* JA     :INOUT    Integer JA(NELT). */
/* A      :INOUT    Real A(NELT). */
/*         These arrays should hold the matrix A in the SLAP Column */
/*         format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* DINV   :OUT      Real DINV(N). */
/*         Upon return this array holds 1./DIAG(A). */

/* *Description */
/*       =================== S L A P Column format ================== */
/*       This routine  requires that  the matrix A  be stored in  the */
/*       SLAP Column format.  In this format the non-zeros are stored */
/*       counting down columns (except  for the diagonal entry, which */
/*       must appear first in each  "column") and  are stored  in the */
/*       real array A.  In other words, for each column in the matrix */
/*       put the diagonal entry in A.  Then put in the other non-zero */
/*       elements going down   the  column (except  the diagonal)  in */
/*       order.  The IA array holds the row  index for each non-zero. */
/*       The JA array holds the offsets into the IA, A arrays for the */
/*       beginning of   each    column.    That  is,    IA(JA(ICOL)), */
/*       A(JA(ICOL)) points to the beginning of the ICOL-th column in */
/*       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the */
/*       end  of   the ICOL-th  column.  Note   that  we  always have */
/*       JA(N+1) = NELT+1, where  N  is the number of columns in  the */
/*       matrix and  NELT   is the number of non-zeros in the matrix. */

/*       Here is an example of the  SLAP Column  storage format for a */
/*       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a */
/*       column): */

/*           5x5 Matrix      SLAP Column format for 5x5 matrix on left. */
/*                              1  2  3    4  5    6  7    8    9 10 11 */
/*       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 */
/*       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 */
/*       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       With the SLAP  format  all  of  the   "inner  loops" of this */
/*       routine should vectorize  on  machines with hardware support */
/*       for vector   gather/scatter  operations.  Your compiler  may */
/*       require a compiler directive to  convince it that  there are */
/*       no  implicit  vector  dependencies.  Compiler directives for */
/*       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are */
/*       supplied with the standard SLAP distribution. */

/* *Precision:           Single Precision */

/* *Cautions: */
/*       This routine assumes that the diagonal of A is all  non-zero */
/*       and that the operation DINV = 1.0/DIAG(A) will not underflow */
/*       or overflow.    This  is done so that the  loop  vectorizes. */
/*       Matricies with zero or near zero or very  large entries will */
/*       have numerical difficulties  and  must  be fixed before this */
/*       routine is called. */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SSDS */

/*         Assume the Diagonal elements are the first in each column. */
/*         This loop should *VECTORIZE*.  If it does not you may have */
/*         to add a compiler directive.  We do not check for a zero */
/*         (or near zero) diagonal element since this would interfere */
/*         with vectorization.  If this makes you nervous put a check */
/*         in!  It will run much slower. */
/* ***FIRST EXECUTABLE STATEMENT  SSDS */
    /* Parameter adjustments */
    --dinv;
    --a;
    --ja;
    --ia;

    /* Function Body */
/* L1: */
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	dinv[icol] = 1.f / a[ja[icol]];
/* L10: */
    }

    return 0;
/* ------------- LAST LINE OF SSDS FOLLOWS ---------------------------- */
} /* ssds_ */

/* DECK SSDSCL */
/* Subroutine */ int ssdscl_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, real *x, real *b, real *dinv, integer *
	job, integer *itol)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j;
    static real di;
    static integer jbgn, jend, icol;

/* ***BEGIN PROLOGUE  SSDSCL */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSDSCL-S), */
/*             SLAP Sparse, Diagonal */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Diagonal Scaling of system Ax = b. */
/*            This routine scales (and unscales) the system Ax = b */
/*            by symmetric diagonal scaling.  The new system is: */
/*             -1/2  -1/2  1/2      -1/2 */
/*            D    AD    (D   x) = D    b */
/*            when scaling is selected with the JOB parameter.  When */
/*            unscaling is selected this process is reversed. */
/*            The true solution is also scaled or unscaled if ITOL is set */
/*            appropriately, see below. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, JOB, ITOL */
/*     REAL    A(NELT), DINV(N) */

/*     CALL SSDSCL( N, NELT, IA, JA, A, ISYM, X, B, DINV, JOB, ITOL ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of elements in arrays IA, JA, and A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Real A(NELT). */
/*         These arrays should hold the matrix A in the SLAP Column */
/*         format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* X      :INOUT    Real X(N). */
/*         Initial guess that will be later used in the iterative */
/*         solution. */
/*         of the scaled system. */
/* B      :INOUT    Real B(N). */
/*         Right hand side vector. */
/* DINV   :OUT      Real DINV(N). */
/*         Upon return this array holds 1./DIAG(A). */
/* JOB    :IN       Integer. */
/*         Flag indicating weather to scale or not.  JOB nonzero means */
/*         do scaling.  JOB = 0 means do unscaling. */
/* ITOL   :IN       Integer. */
/*         Flag indicating what type of error estimation to do in the */
/*         iterative method.  When ITOL = 11 the exact solution from */
/*         common block solblk will be used.  When the system is scaled */
/*         then the true solution must also be scaled.  If ITOL is not */
/*         11 then this vector is not referenced. */

/* *Common Blocks: */
/* SOLN    :INOUT   Real SOLN(N).  COMMON BLOCK /SOLBLK/ */
/*         The true solution, SOLN, is scaled (or unscaled) if ITOL is */
/*         set to 11, see above. */

/* *Description */
/*       =================== S L A P Column format ================== */
/*       This routine  requires that  the matrix A  be stored in  the */
/*       SLAP Column format.  In this format the non-zeros are stored */
/*       counting down columns (except  for the diagonal entry, which */
/*       must appear first in each  "column") and  are stored  in the */
/*       real array A.  In other words, for each column in the matrix */
/*       put the diagonal entry in A.  Then put in the other non-zero */
/*       elements going down   the  column (except  the diagonal)  in */
/*       order.  The IA array holds the row  index for each non-zero. */
/*       The JA array holds the offsets into the IA, A arrays for the */
/*       beginning of   each    column.    That  is,    IA(JA(ICOL)), */
/*       A(JA(ICOL)) points to the beginning of the ICOL-th column in */
/*       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the */
/*       end  of   the ICOL-th  column.  Note   that  we  always have */
/*       JA(N+1) = NELT+1, where  N  is the number of columns in  the */
/*       matrix and  NELT   is the number of non-zeros in the matrix. */

/*       Here is an example of the  SLAP Column  storage format for a */
/*       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a */
/*       column): */

/*           5x5 Matrix      SLAP Column format for 5x5 matrix on left. */
/*                              1  2  3    4  5    6  7    8    9 10 11 */
/*       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 */
/*       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 */
/*       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       With the SLAP  format  all  of  the   "inner  loops" of this */
/*       routine should vectorize  on  machines with hardware support */
/*       for vector   gather/scatter  operations.  Your compiler  may */
/*       require a compiler directive to  convince it that  there are */
/*       no  implicit  vector  dependencies.  Compiler directives for */
/*       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are */
/*       supplied with the standard SLAP distribution. */

/* *Precision:           Single Precision */

/* *Cautions: */
/*       This routine assumes that the diagonal of A is all  non-zero */
/*       and that the operation DINV = 1.0/DIAG(A)  will  not  under- */
/*       flow or overflow. This is done so that the loop  vectorizes. */
/*       Matricies with zero or near zero or very  large entries will */
/*       have numerical difficulties  and  must  be fixed before this */
/*       routine is called. */

/* *See Also: */
/*       SSDCG */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SSDSCL */

/*         SCALING... */

    /* Parameter adjustments */
    --dinv;
    --b;
    --x;
    --a;
    --ja;
    --ia;

    /* Function Body */
    if (*job != 0) {
	i__1 = *n;
	for (icol = 1; icol <= i__1; ++icol) {
	    dinv[icol] = 1.f / sqrt(a[ja[icol]]);
/* L10: */
	}
    } else {

/*         UNSCALING... */

	i__1 = *n;
	for (icol = 1; icol <= i__1; ++icol) {
	    dinv[icol] = 1.f / dinv[icol];
/* L15: */
	}
    }

    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	jbgn = ja[icol];
	jend = ja[icol + 1] - 1;
	di = dinv[icol];
	i__2 = jend;
	for (j = jbgn; j <= i__2; ++j) {
	    a[j] = dinv[ia[j]] * a[j] * di;
/* L20: */
	}
/* L30: */
    }

    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	b[icol] *= dinv[icol];
	x[icol] /= dinv[icol];
/* L40: */
    }

/*         Check to see if we need to scale the "true solution" as well. */

    if (*itol == 11) {
	i__1 = *n;
	for (icol = 1; icol <= i__1; ++icol) {
	    solblk_1.soln[icol - 1] /= dinv[icol];
/* L50: */
	}
    }

    return 0;
} /* ssdscl_ */

/* DECK SSD2S */
/* Subroutine */ int ssd2s_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, real *dinv)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, k, kbgn, kend;

/* ***BEGIN PROLOGUE  SSD2S */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSD2S-S), */
/*             SLAP Sparse, Diagonal */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Diagonal Scaling Preconditioner SLAP Normal Eqns Set Up. */
/*            Routine to compute the inverse of the diagonal of the */
/*            matrix A*A'.  Where A is stored in SLAP-Column format. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM */
/*     REAL    A(NELT), DINV(N) */

/*     CALL SSD2S( N, NELT, IA, JA, A, ISYM, DINV ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of elements in arrays IA, JA, and A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Real A(NELT). */
/*         These arrays should hold the matrix A in the SLAP Column */
/*         format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* DINV   :OUT      Real DINV(N). */
/*         Upon return this array holds 1./DIAG(A*A'). */

/* *Description */
/*       =================== S L A P Column format ================== */
/*       This routine  requires that  the matrix A  be stored in  the */
/*       SLAP Column format.  In this format the non-zeros are stored */
/*       counting down columns (except  for the diagonal entry, which */
/*       must appear first in each  "column") and  are stored  in the */
/*       real array A.  In other words, for each column in the matrix */
/*       put the diagonal entry in A.  Then put in the other non-zero */
/*       elements going down   the  column (except  the diagonal)  in */
/*       order.  The IA array holds the row  index for each non-zero. */
/*       The JA array holds the offsets into the IA, A arrays for the */
/*       beginning of   each    column.    That  is,    IA(JA(ICOL)), */
/*       A(JA(ICOL)) points to the beginning of the ICOL-th column in */
/*       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the */
/*       end  of   the ICOL-th  column.  Note   that  we  always have */
/*       JA(N+1) = NELT+1, where  N  is the number of columns in  the */
/*       matrix and  NELT   is the number of non-zeros in the matrix. */

/*       Here is an example of the  SLAP Column  storage format for a */
/*       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a */
/*       column): */

/*           5x5 Matrix      SLAP Column format for 5x5 matrix on left. */
/*                              1  2  3    4  5    6  7    8    9 10 11 */
/*       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 */
/*       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 */
/*       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       With the SLAP  format  all  of  the   "inner  loops" of this */
/*       routine should vectorize  on  machines with hardware support */
/*       for vector   gather/scatter  operations.  Your compiler  may */
/*       require a compiler directive to  convince it that  there are */
/*       no  implicit  vector  dependencies.  Compiler directives for */
/*       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are */
/*       supplied with the standard SLAP distribution. */

/* *Precision:           Single Precision */

/* *Cautions: */
/*       This routine assumes that the diagonal of A is all  non-zero */
/*       and that the operation DINV = 1.0/DIAG(A*A') will not under- */
/*       flow or overflow. This is done so that the loop  vectorizes. */
/*       Matricies with zero or near zero or very  large entries will */
/*       have numerical difficulties  and  must  be fixed before this */
/*       routine is called. */

/* *See Also: */
/*       SSDCGN */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SSD2S */

/* ***FIRST EXECUTABLE STATEMENT  SSD2S */
    /* Parameter adjustments */
    --dinv;
    --a;
    --ja;
    --ia;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dinv[i__] = 0.f;
/* L10: */
    }

/*         Loop over each column. */
/* VD$R NOCONCUR */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	kbgn = ja[i__];
	kend = ja[i__ + 1] - 1;

/*         Add in the contributions for each row that has a non-zero */
/*         in this column. */
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
	i__2 = kend;
	for (k = kbgn; k <= i__2; ++k) {
/* Computing 2nd power */
	    r__1 = a[k];
	    dinv[ia[k]] += r__1 * r__1;
/* L20: */
	}
	if (*isym == 1) {

/*         Lower triangle stored by columns => upper triangle stored by */
/*         rows with Diagonal being the first entry.  Loop across the */
/*         rest of the row. */
	    ++kbgn;
	    if (kbgn <= kend) {
		i__2 = kend;
		for (k = kbgn; k <= i__2; ++k) {
/* Computing 2nd power */
		    r__1 = a[k];
		    dinv[i__] += r__1 * r__1;
/* L30: */
		}
	    }
	}
/* L40: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dinv[i__] = 1.f / dinv[i__];
/* L50: */
    }

    return 0;
/* ------------- LAST LINE OF SSD2S FOLLOWS ---------------------------- */
} /* ssd2s_ */

/* DECK SS2LT */
/* Subroutine */ int ss2lt_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, integer *nel, integer *iel, integer *jel, 
	real *el)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, jbgn, jend, icol;

/* ***BEGIN PROLOGUE  SS2LT */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SS2LT-S), */
/*             Linear system, SLAP Sparse, Lower Triangle */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Lower Triangle Preconditioner SLAP Set Up. */
/*            Routine to store the lower triangle of a matrix stored */
/*            in the Slap Column format. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM */
/*     INTEGER NEL, IEL(N+1), JEL(NEL), NROW(N) */
/*     REAL    A(NELT), EL(NEL) */

/*     CALL SS2LT( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of non-zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Real A(NELT). */
/*         These arrays should hold the matrix A in the SLAP Column */
/*         format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* NEL    :OUT      Integer. */
/*         Number of non-zeros in the lower triangle of A.   Also */
/*         coresponds to the length of the JEL, EL arrays. */
/* IEL    :OUT      Integer IEL(N+1). */
/* JEL    :OUT      Integer JEL(NEL). */
/* EL     :OUT      Real     EL(NEL). */
/*         IEL, JEL, EL contain the lower triangle of the A matrix */
/*         stored in SLAP Column format.  See "Description", below */
/*         for more details bout the SLAP Column format. */

/* *Description */
/*       =================== S L A P Column format ================== */
/*       This routine  requires that  the matrix A  be stored in  the */
/*       SLAP Column format.  In this format the non-zeros are stored */
/*       counting down columns (except  for the diagonal entry, which */
/*       must appear first in each  "column") and  are stored  in the */
/*       real array A.  In other words, for each column in the matrix */
/*       put the diagonal entry in A.  Then put in the other non-zero */
/*       elements going down   the  column (except  the diagonal)  in */
/*       order.  The IA array holds the row  index for each non-zero. */
/*       The JA array holds the offsets into the IA, A arrays for the */
/*       beginning of   each    column.    That  is,    IA(JA(ICOL)), */
/*       A(JA(ICOL)) points to the beginning of the ICOL-th column in */
/*       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the */
/*       end  of   the ICOL-th  column.  Note   that  we  always have */
/*       JA(N+1) = NELT+1, where  N  is the number of columns in  the */
/*       matrix and  NELT   is the number of non-zeros in the matrix. */

/*       Here is an example of the  SLAP Column  storage format for a */
/*       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a */
/*       column): */

/*           5x5 Matrix      SLAP Column format for 5x5 matrix on left. */
/*                              1  2  3    4  5    6  7    8    9 10 11 */
/*       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 */
/*       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 */
/*       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/* *Precision:           Single Precision */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SS2LT */
/* ***FIRST EXECUTABLE STATEMENT  SS2LT */
    /* Parameter adjustments */
    --el;
    --a;
    --ja;
    --ia;
    --jel;
    --iel;

    /* Function Body */
    if (*isym == 0) {

/*         The matrix is stored non-symmetricly.  Pick out the lower */
/*         triangle. */

	*nel = 0;
	i__1 = *n;
	for (icol = 1; icol <= i__1; ++icol) {
	    jel[icol] = *nel + 1;
	    jbgn = ja[icol];
	    jend = ja[icol + 1] - 1;
/* VD$ NOVECTOR */
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		if (ia[j] >= icol) {
		    ++(*nel);
		    iel[*nel] = ia[j];
		    el[*nel] = a[j];
		}
/* L10: */
	    }
/* L20: */
	}
	jel[*n + 1] = *nel + 1;
    } else {

/*         The matrix is symmetric and only the lower triangle is */
/*         stored.  Copy it to IEL, JEL, EL. */

	*nel = *nelt;
	i__1 = *nelt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    iel[i__] = ia[i__];
	    el[i__] = a[i__];
/* L30: */
	}
	i__1 = *n + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    jel[i__] = ja[i__];
/* L40: */
	}
    }
    return 0;
/* ------------- LAST LINE OF SS2LT FOLLOWS ---------------------------- */
} /* ss2lt_ */

/* DECK SSICS */
/* Subroutine */ int ssics_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, integer *nel, integer *iel, integer *jel, 
	real *el, real *d__, real *r__, integer *iwarn)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Local variables */
    static integer i__, j, ic, ir, irr, ibgn, jbgn, jend, iend, icol, irow, 
	    icbgn, icend, irbgn, irend;
    static real eltmp;
    static integer jeltmp;
    extern /* Subroutine */ int xerrwv_(char *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, real *, real *, 
	    ftnlen);

/* ***BEGIN PROLOGUE  SSICS */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSICS-S), */
/*             Linear system, SLAP Sparse, Iterative Precondition */
/*             Incomplete Cholesky Factorization. */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Incompl Cholesky Decomposition Preconditioner SLAP Set Up. */
/*            Routine to generate the Incomplete Cholesky decomposition, */
/*            L*D*L-trans, of  a symmetric positive definite  matrix, A, */
/*            which  is stored  in  SLAP Column format.  The  unit lower */
/*            triangular matrix L is  stored by rows, and the inverse of */
/*            the diagonal matrix D is stored. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM */
/*     INTEGER NEL, IEL(NEL), JEL(N+1), IWARN */
/*     REAL    A(NELT), EL(NEL), D(N), R(N) */

/*     CALL SSICS( N, NELT, IA, JA, A, ISYM, NEL, IEL, JEL, EL, D, R, */
/*    $    IWARN ) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of elements in arrays IA, JA, and A. */
/* IA     :INOUT    Integer IA(NELT). */
/* JA     :INOUT    Integer JA(NELT). */
/* A      :INOUT    Real A(NELT). */
/*         These arrays should hold the matrix A in the SLAP Column */
/*         format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* NEL    :OUT      Integer. */
/*         Number of non-zeros in the lower triangle of A.   Also */
/*         coresponds to the length of the JEL, EL arrays. */
/* IEL    :OUT      Integer IEL(N+1). */
/* JEL    :OUT      Integer JEL(NEL). */
/* EL     :OUT      Real     EL(NEL). */
/*         IEL, JEL, EL contain the unit lower triangular factor  of the */
/*         incomplete decomposition   of the A  matrix  stored  in  SLAP */
/*         Row format.   The Diagonal of   ones   *IS*   stored.     See */
/*         "Description", below for more details about the SLAP Row fmt. */
/* D      :OUT      Real D(N) */
/*         Upon return this array holds D(I) = 1./DIAG(A). */
/* R      :WORK     Real R(N). */
/*         Temporary real workspace needed for the factorization. */
/* IWARN  :OUT      Integer. */
/*         This is a warning variable and is zero if the IC factoriza- */
/*         tion goes well.  It is set to the row index corresponding to */
/*         the last zero pivot found.  See "Description", below. */

/* *Description */
/*       =================== S L A P Column format ================== */
/*       This routine  requires that  the matrix A  be stored in  the */
/*       SLAP Column format.  In this format the non-zeros are stored */
/*       counting down columns (except  for the diagonal entry, which */
/*       must appear first in each  "column") and  are stored  in the */
/*       real array A.  In other words, for each column in the matrix */
/*       put the diagonal entry in A.  Then put in the other non-zero */
/*       elements going down   the  column (except  the diagonal)  in */
/*       order.  The IA array holds the row  index for each non-zero. */
/*       The JA array holds the offsets into the IA, A arrays for the */
/*       beginning of   each    column.    That  is,    IA(JA(ICOL)), */
/*       A(JA(ICOL)) points to the beginning of the ICOL-th column in */
/*       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the */
/*       end  of   the ICOL-th  column.  Note   that  we  always have */
/*       JA(N+1) = NELT+1, where  N  is the number of columns in  the */
/*       matrix and  NELT   is the number of non-zeros in the matrix. */

/*       Here is an example of the  SLAP Column  storage format for a */
/*       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a */
/*       column): */

/*           5x5 Matrix      SLAP Column format for 5x5 matrix on left. */
/*                              1  2  3    4  5    6  7    8    9 10 11 */
/*       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 */
/*       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 */
/*       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       ==================== S L A P Row format ==================== */
/*       This routine requires  that the matrix A  be  stored  in the */
/*       SLAP  Row format.   In this format  the non-zeros are stored */
/*       counting across  rows (except for the diagonal  entry, which */
/*       must appear first in each "row") and  are stored in the real */
/*       array A.  In other words, for each row in the matrix put the */
/*       diagonal entry in  A.   Then   put  in the   other  non-zero */
/*       elements   going  across the  row (except   the diagonal) in */
/*       order.   The  JA array  holds   the column   index for  each */
/*       non-zero.   The IA  array holds the  offsets into  the JA, A */
/*       arrays  for   the   beginning  of   each  row.   That    is, */
/*       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the */
/*       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1) */
/*       points to the  end of the  IROW-th row.  Note that we always */
/*       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in */
/*       the matrix  and NELT  is the  number   of  non-zeros in  the */
/*       matrix. */

/*       Here is an example of the SLAP Row storage format for a  5x5 */
/*       Matrix (in the A and JA arrays '|' denotes the end of a row): */

/*           5x5 Matrix         SLAP Row format for 5x5 matrix on left. */
/*                              1  2  3    4  5    6  7    8    9 10 11 */
/*       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53 */
/*       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 */
/*       | 0  0 33  0 35|  IA:  1  4  6    8  9   12 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       With the SLAP  format some  of  the   "inner  loops" of this */
/*       routine should vectorize  on  machines with hardware support */
/*       for vector   gather/scatter  operations.  Your compiler  may */
/*       require a compiler directive to  convince it that  there are */
/*       no  implicit  vector  dependencies.  Compiler directives for */
/*       the Alliant    FX/Fortran and CRI   CFT/CFT77 compilers  are */
/*       supplied with the standard SLAP distribution. */

/*       The IC  factorization is not  alway exist for SPD matricies. */
/*       In the event that a zero pivot is found it is set  to be 1.0 */
/*       and the factorization procedes.   The integer variable IWARN */
/*       is set to the last row where the Diagonal was fudged.  This */
/*       eventuality hardly ever occurs in practice */

/* *Precision:           Single Precision */

/* *See Also: */
/*       SCG, SSICCG */
/* ***REFERENCES  1. Gene Golub & Charles Van Loan, "Matrix Computations", */
/*                 John Hopkins University Press; 3 (1983) IBSN */
/*                 0-8018-3010-9. */
/* ***ROUTINES CALLED  XERRWV */
/* ***END PROLOGUE  SSICS */

/*         Set the lower triangle in IEL, JEL, EL */
/* ***FIRST EXECUTABLE STATEMENT  SSICS */
    /* Parameter adjustments */
    --r__;
    --d__;
    --a;
    --ja;
    --ia;
    --el;
    --jel;
    --iel;

    /* Function Body */
    *iwarn = 0;

/*         All matrix elements stored in IA, JA, A.  Pick out the lower */
/*         triangle (making sure that the Diagonal of EL is one) and */
/*         store by rows. */

    *nel = 1;
    iel[1] = 1;
    jel[1] = 1;
    el[1] = 1.f;
    d__[1] = a[1];
/* VD$R NOCONCUR */
    i__1 = *n;
    for (irow = 2; irow <= i__1; ++irow) {
/*         Put in the Diagonal. */
	++(*nel);
	iel[irow] = *nel;
	jel[*nel] = irow;
	el[*nel] = 1.f;
	d__[irow] = a[ja[irow]];

/*         Look in all the lower triangle columns for a matching row. */
/*         Since the matrix is symmetric, we can look across the */
/*         irow-th row by looking down the irow-th column (if it is */
/*         stored ISYM=0)... */
	if (*isym == 0) {
	    icbgn = ja[irow];
	    icend = ja[irow + 1] - 1;
	} else {
	    icbgn = 1;
	    icend = irow - 1;
	}
	i__2 = icend;
	for (ic = icbgn; ic <= i__2; ++ic) {
	    if (*isym == 0) {
		icol = ia[ic];
		if (icol >= irow) {
		    goto L20;
		}
	    } else {
		icol = ic;
	    }
	    jbgn = ja[icol] + 1;
	    jend = ja[icol + 1] - 1;
	    if (jbgn <= jend && ia[jend] >= irow) {
/* VD$ NOVECTOR */
		i__3 = jend;
		for (j = jbgn; j <= i__3; ++j) {
		    if (ia[j] == irow) {
			++(*nel);
			jel[*nel] = icol;
			el[*nel] = a[j];
			goto L20;
		    }
/* L10: */
		}
	    }
L20:
	    ;
	}
/* L30: */
    }
    iel[*n + 1] = *nel + 1;

/*         Sort ROWS of lower triangle into descending order (count out */
/*         along rows out from Diagonal). */

    i__1 = *n;
    for (irow = 2; irow <= i__1; ++irow) {
	ibgn = iel[irow] + 1;
	iend = iel[irow + 1] - 1;
	if (ibgn < iend) {
	    i__2 = iend - 1;
	    for (i__ = ibgn; i__ <= i__2; ++i__) {
/* VD$ NOVECTOR */
		i__3 = iend;
		for (j = i__ + 1; j <= i__3; ++j) {
		    if (jel[i__] > jel[j]) {
			jeltmp = jel[j];
			jel[j] = jel[i__];
			jel[i__] = jeltmp;
			eltmp = el[j];
			el[j] = el[i__];
			el[i__] = eltmp;
		    }
/* L40: */
		}
/* L50: */
	    }
	}
/* L60: */
    }

/*         Perform the Incomplete Cholesky decomposition by looping */
/*         over the rows. */
/*         Scale the first column.  Use the structure of A to pick out */
/*         the rows with something in column 1. */

    irbgn = ja[1] + 1;
    irend = ja[2] - 1;
    i__1 = irend;
    for (irr = irbgn; irr <= i__1; ++irr) {
	ir = ia[irr];
/*         Find the index into EL for EL(1,IR). */
/*         Hint: it's the second entry. */
	i__ = iel[ir] + 1;
	el[i__] /= d__[1];
/* L65: */
    }

    i__1 = *n;
    for (irow = 2; irow <= i__1; ++irow) {

/*         Update the IROW-th diagonal. */

	i__2 = irow - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__] = 0.f;
/* L66: */
	}
	ibgn = iel[irow] + 1;
	iend = iel[irow + 1] - 1;
	if (ibgn <= iend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
	    i__2 = iend;
	    for (i__ = ibgn; i__ <= i__2; ++i__) {
		r__[jel[i__]] = el[i__] * d__[jel[i__]];
		d__[irow] -= el[i__] * r__[jel[i__]];
/* L70: */
	    }

/*         Check to see if we gota problem with the diagonal. */

	    if (d__[irow] <= 0.f) {
		if (*iwarn == 0) {
		    *iwarn = irow;
		}
		d__[irow] = 1.f;
	    }
	}

/*         Update each EL(IROW+1:N,IROW), if there are any. */
/*         Use the structure of A to determine the Non-zero elements */
/*         of the IROW-th column of EL. */

	irbgn = ja[irow];
	irend = ja[irow + 1] - 1;
	i__2 = irend;
	for (irr = irbgn; irr <= i__2; ++irr) {
	    ir = ia[irr];
	    if (ir <= irow) {
		goto L100;
	    }
/*         Find the index into EL for EL(IR,IROW) */
	    ibgn = iel[ir] + 1;
	    iend = iel[ir + 1] - 1;
	    if (jel[ibgn] > irow) {
		goto L100;
	    }
	    i__3 = iend;
	    for (i__ = ibgn; i__ <= i__3; ++i__) {
		if (jel[i__] == irow) {
		    icend = iend;
L91:
		    if (jel[icend] >= irow) {
			--icend;
			goto L91;
		    }
/*         Sum up the EL(IR,1:IROW-1)*R(1:IROW-1) contributions. */
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
		    i__4 = icend;
		    for (ic = ibgn; ic <= i__4; ++ic) {
			el[i__] -= el[ic] * r__[jel[ic]];
/* L80: */
		    }
		    el[i__] /= d__[irow];
		    goto L100;
		}
/* L90: */
	    }

/*         If we get here, we have real problems... */
	    xerrwv_("SSICS -- A and EL data structure mismatch in row (i1)", &
		    c__53, &c__1, &c__2, &c__1, &irow, &c__0, &c__0, &c_b44, &
		    c_b44, (ftnlen)53);
L100:
	    ;
	}
/* L110: */
    }

/*         Replace diagonals by their inverses. */

/* VD$ CONCUR */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = 1.f / d__[i__];
/* L120: */
    }
    return 0;
/* ------------- LAST LINE OF SSICS FOLLOWS ---------------------------- */
} /* ssics_ */

/* DECK SSILUS */
/* Subroutine */ int ssilus_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, integer *nl, integer *il, integer *jl, 
	real *l, real *dinv, integer *nu, integer *iu, integer *ju, real *u, 
	integer *nrow, integer *ncol)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, kc, kr, ibgn, jbgn, jend, iend, icol, indx;
    static real temp;
    static integer irow, indx1, indx2, itemp, jtemp, indxc1, indxc2, indxr1, 
	    indxr2;

/* ***BEGIN PROLOGUE  SSILUS */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSILUS-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Incomplete LU Factorization */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Incomplete LU Decomposition Preconditioner SLAP Set Up. */
/*            Routine to generate the incomplete LDU decomposition of a */
/*            matrix.  The  unit lower triangular factor L is stored by */
/*            rows and the  unit upper triangular factor U is stored by */
/*            columns.  The inverse of the diagonal matrix D is stored. */
/*            No fill in is allowed. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM */
/*     INTEGER NL, IL(N+1), JL(NL), NU, IU(N+1), JU(NU) */
/*     INTEGER NROW(N), NCOL(N) */
/*     REAL    A(NELT), L(NL), U(NU), DINV(N) */

/*     CALL SSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, */
/*    $    DINV, NU, IU, JU, U, NROW, NCOL ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of elements in arrays IA, JA, and A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Real A(NELT). */
/*         These arrays should hold the matrix A in the SLAP Column */
/*         format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* NL     :OUT      Integer. */
/*         Number of non-zeros in the EL array. */
/* IL     :OUT      Integer IL(N+1). */
/* JL     :OUT      Integer JL(NL). */
/* L      :OUT      Real     L(NL). */
/*         IL, JL, L  contain the unit ower  triangular factor of  the */
/*         incomplete decomposition  of some  matrix stored  in   SLAP */
/*         Row format.     The   Diagonal  of ones  *IS*  stored.  See */
/*         "DESCRIPTION", below for more details about the SLAP format. */
/* NU     :OUT      Integer. */
/*         Number of non-zeros in the U array. */
/* IU     :OUT      Integer IU(N+1). */
/* JU     :OUT      Integer JU(NU). */
/* U      :OUT      Real     U(NU). */
/*         IU, JU, U contain   the unit upper triangular factor of the */
/*         incomplete  decomposition    of some matrix  stored in SLAP */
/*         Column  format.   The Diagonal of ones   *IS*  stored.  See */
/*         "Description", below  for  more  details  about  the   SLAP */
/*         format. */
/* NROW   :WORK     Integer NROW(N). */
/*         NROW(I) is the number of non-zero elements in the I-th row */
/*         of L. */
/* NCOL   :WORK     Integer NCOL(N). */
/*         NCOL(I) is the number of non-zero elements in the I-th */
/*         column of U. */

/* *Description */
/*       IL, JL, L should contain the unit  lower triangular factor of */
/*       the incomplete decomposition of the A matrix  stored in SLAP */
/*       Row format.  IU, JU, U should contain  the unit upper factor */
/*       of the  incomplete decomposition of  the A matrix  stored in */
/*       SLAP Column format This ILU factorization can be computed by */
/*       the SSILUS routine.  The diagonals (which is all one's) are */
/*       stored. */

/*       =================== S L A P Column format ================== */
/*       This routine  requires that  the matrix A  be stored in  the */
/*       SLAP Column format.  In this format the non-zeros are stored */
/*       counting down columns (except  for the diagonal entry, which */
/*       must appear first in each  "column") and  are stored  in the */
/*       real array A.  In other words, for each column in the matrix */
/*       put the diagonal entry in A.  Then put in the other non-zero */
/*       elements going down   the  column (except  the diagonal)  in */
/*       order.  The IA array holds the row  index for each non-zero. */
/*       The JA array holds the offsets into the IA, A arrays for the */
/*       beginning of   each    column.    That  is,    IA(JA(ICOL)), */
/*       A(JA(ICOL)) points to the beginning of the ICOL-th column in */
/*       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the */
/*       end  of   the ICOL-th  column.  Note   that  we  always have */
/*       JA(N+1) = NELT+1, where  N  is the number of columns in  the */
/*       matrix and  NELT   is the number of non-zeros in the matrix. */

/*       Here is an example of the  SLAP Column  storage format for a */
/*       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a */
/*       column): */

/*           5x5 Matrix      SLAP Column format for 5x5 matrix on left. */
/*                              1  2  3    4  5    6  7    8    9 10 11 */
/*       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35 */
/*       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 */
/*       | 0  0 33  0 35|  JA:  1  4  6    8  9   12 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       ==================== S L A P Row format ==================== */
/*       This routine requires  that the matrix A  be  stored  in the */
/*       SLAP  Row format.   In this format  the non-zeros are stored */
/*       counting across  rows (except for the diagonal  entry, which */
/*       must appear first in each "row") and  are stored in the real */
/*       array A.  In other words, for each row in the matrix put the */
/*       diagonal entry in  A.   Then   put  in the   other  non-zero */
/*       elements   going  across the  row (except   the diagonal) in */
/*       order.   The  JA array  holds   the column   index for  each */
/*       non-zero.   The IA  array holds the  offsets into  the JA, A */
/*       arrays  for   the   beginning  of   each  row.   That    is, */
/*       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the */
/*       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1) */
/*       points to the  end of the  IROW-th row.  Note that we always */
/*       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in */
/*       the matrix  and NELT  is the  number   of  non-zeros in  the */
/*       matrix. */

/*       Here is an example of the SLAP Row storage format for a  5x5 */
/*       Matrix (in the A and JA arrays '|' denotes the end of a row): */

/*           5x5 Matrix         SLAP Row format for 5x5 matrix on left. */
/*                              1  2  3    4  5    6  7    8    9 10 11 */
/*       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53 */
/*       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3 */
/*       | 0  0 33  0 35|  IA:  1  4  6    8  9   12 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/* *Precision:           Single Precision */
/* *See Also: */
/*       SILUR */
/* ***REFERENCES  1. Gene Golub & Charles Van Loan, "Matrix Computations", */
/*                 John Hopkins University Press; 3 (1983) IBSN */
/*                 0-8018-3010-9. */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SSILUS */

/*         Count number of elements in each row of the lower triangle. */
/* ***FIRST EXECUTABLE STATEMENT  SSILUS */
    /* Parameter adjustments */
    --ncol;
    --nrow;
    --dinv;
    --a;
    --ja;
    --ia;
    --l;
    --jl;
    --il;
    --u;
    --ju;
    --iu;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	nrow[i__] = 0;
	ncol[i__] = 0;
/* L10: */
    }
/* VD$R NOCONCUR */
/* VD$R NOVECTOR */
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	jbgn = ja[icol] + 1;
	jend = ja[icol + 1] - 1;
	if (jbgn <= jend) {
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		if (ia[j] < icol) {
		    ++ncol[icol];
		} else {
		    ++nrow[ia[j]];
		    if (*isym != 0) {
			++ncol[ia[j]];
		    }
		}
/* L20: */
	    }
	}
/* L30: */
    }
    ju[1] = 1;
    il[1] = 1;
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	il[icol + 1] = il[icol] + nrow[icol];
	ju[icol + 1] = ju[icol] + ncol[icol];
	nrow[icol] = il[icol];
	ncol[icol] = ju[icol];
/* L40: */
    }

/*         Copy the matrix A into the L and U structures. */
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	dinv[icol] = a[ja[icol]];
	jbgn = ja[icol] + 1;
	jend = ja[icol + 1] - 1;
	if (jbgn <= jend) {
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		irow = ia[j];
		if (irow < icol) {
/*         Part of the upper triangle. */
		    iu[ncol[icol]] = irow;
		    u[ncol[icol]] = a[j];
		    ++ncol[icol];
		} else {
/*         Part of the lower triangle (stored by row). */
		    jl[nrow[irow]] = icol;
		    l[nrow[irow]] = a[j];
		    ++nrow[irow];
		    if (*isym != 0) {
/*         Symmetric...Copy lower triangle into upper triangle as well. */
			iu[ncol[irow]] = icol;
			u[ncol[irow]] = a[j];
			++ncol[irow];
		    }
		}
/* L50: */
	    }
	}
/* L60: */
    }

/*         Sort the rows of L and the columns of U. */
    i__1 = *n;
    for (k = 2; k <= i__1; ++k) {
	jbgn = ju[k];
	jend = ju[k + 1] - 1;
	if (jbgn < jend) {
	    i__2 = jend - 1;
	    for (j = jbgn; j <= i__2; ++j) {
		i__3 = jend;
		for (i__ = j + 1; i__ <= i__3; ++i__) {
		    if (iu[j] > iu[i__]) {
			itemp = iu[j];
			iu[j] = iu[i__];
			iu[i__] = itemp;
			temp = u[j];
			u[j] = u[i__];
			u[i__] = temp;
		    }
/* L70: */
		}
/* L80: */
	    }
	}
	ibgn = il[k];
	iend = il[k + 1] - 1;
	if (ibgn < iend) {
	    i__2 = iend - 1;
	    for (i__ = ibgn; i__ <= i__2; ++i__) {
		i__3 = iend;
		for (j = i__ + 1; j <= i__3; ++j) {
		    if (jl[i__] > jl[j]) {
			jtemp = ju[i__];
			ju[i__] = ju[j];
			ju[j] = jtemp;
			temp = l[i__];
			l[i__] = l[j];
			l[j] = temp;
		    }
/* L90: */
		}
/* L100: */
	    }
	}
/* L110: */
    }

/*         Perform the incomplete LDU decomposition. */
    i__1 = *n;
    for (i__ = 2; i__ <= i__1; ++i__) {

/*           I-th row of L */
	indx1 = il[i__];
	indx2 = il[i__ + 1] - 1;
	if (indx1 > indx2) {
	    goto L200;
	}
	i__2 = indx2;
	for (indx = indx1; indx <= i__2; ++indx) {
	    if (indx == indx1) {
		goto L180;
	    }
	    indxr1 = indx1;
	    indxr2 = indx - 1;
	    indxc1 = ju[jl[indx]];
	    indxc2 = ju[jl[indx] + 1] - 1;
	    if (indxc1 > indxc2) {
		goto L180;
	    }
L160:
	    kr = jl[indxr1];
L170:
	    kc = iu[indxc1];
	    if (kr > kc) {
		++indxc1;
		if (indxc1 <= indxc2) {
		    goto L170;
		}
	    } else if (kr < kc) {
		++indxr1;
		if (indxr1 <= indxr2) {
		    goto L160;
		}
	    } else if (kr == kc) {
		l[indx] -= l[indxr1] * dinv[kc] * u[indxc1];
		++indxr1;
		++indxc1;
		if (indxr1 <= indxr2 && indxc1 <= indxc2) {
		    goto L160;
		}
	    }
L180:
	    l[indx] /= dinv[jl[indx]];
/* L190: */
	}

/*         ith column of u */
L200:
	indx1 = ju[i__];
	indx2 = ju[i__ + 1] - 1;
	if (indx1 > indx2) {
	    goto L260;
	}
	i__2 = indx2;
	for (indx = indx1; indx <= i__2; ++indx) {
	    if (indx == indx1) {
		goto L240;
	    }
	    indxc1 = indx1;
	    indxc2 = indx - 1;
	    indxr1 = il[iu[indx]];
	    indxr2 = il[iu[indx] + 1] - 1;
	    if (indxr1 > indxr2) {
		goto L240;
	    }
L210:
	    kr = jl[indxr1];
L220:
	    kc = iu[indxc1];
	    if (kr > kc) {
		++indxc1;
		if (indxc1 <= indxc2) {
		    goto L220;
		}
	    } else if (kr < kc) {
		++indxr1;
		if (indxr1 <= indxr2) {
		    goto L210;
		}
	    } else if (kr == kc) {
		u[indx] -= l[indxr1] * dinv[kc] * u[indxc1];
		++indxr1;
		++indxc1;
		if (indxr1 <= indxr2 && indxc1 <= indxc2) {
		    goto L210;
		}
	    }
L240:
	    u[indx] /= dinv[iu[indx]];
/* L250: */
	}

/*         ith diagonal element */
L260:
	indxr1 = il[i__];
	indxr2 = il[i__ + 1] - 1;
	if (indxr1 > indxr2) {
	    goto L300;
	}
	indxc1 = ju[i__];
	indxc2 = ju[i__ + 1] - 1;
	if (indxc1 > indxc2) {
	    goto L300;
	}
L270:
	kr = jl[indxr1];
L280:
	kc = iu[indxc1];
	if (kr > kc) {
	    ++indxc1;
	    if (indxc1 <= indxc2) {
		goto L280;
	    }
	} else if (kr < kc) {
	    ++indxr1;
	    if (indxr1 <= indxr2) {
		goto L270;
	    }
	} else if (kr == kc) {
	    dinv[i__] -= l[indxr1] * dinv[kc] * u[indxc1];
	    ++indxr1;
	    ++indxc1;
	    if (indxr1 <= indxr2 && indxc1 <= indxc2) {
		goto L270;
	    }
	}

L300:
	;
    }

/*         replace diagonal lts by their inverses. */
/* VD$ VECTOR */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dinv[i__] = 1.f / dinv[i__];
/* L430: */
    }

    return 0;
/* ------------- LAST LINE OF SSILUS FOLLOWS ---------------------------- */
} /* ssilus_ */

