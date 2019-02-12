/* smvops.f -- translated by f2c (version 20100827).
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

/* DECK SSMV */
/* Subroutine */ int ssmv_(integer *n, real *x, real *y, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, ibgn, iend, jbgn, jend, icol, irow;

/* ***BEGIN PROLOGUE  SSMV */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSMV-S), */
/*             Matrix Vector Multiply, Sparse */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP Column Format Sparse Matrix Vector Product. */
/*            Routine to calculate the sparse matrix vector product: */
/*            Y = A*X. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER  N, NELT, IA(NELT), JA(N+1), ISYM */
/*     REAL     X(N), Y(N), A(NELT) */

/*     CALL SSMV(N, X, Y, NELT, IA, JA, A, ISYM ) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* X      :IN       Real X(N). */
/*         The vector that should be multiplied by the matrix. */
/* Y      :OUT      Real Y(N). */
/*         The product of the matrix and the vector. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(N+1). */
/* A      :IN       Integer A(NELT). */
/*         These arrays should hold the matrix A in the SLAP Column */
/*         format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */

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

/*       With  the SLAP  format  the "inner  loops" of  this  routine */
/*       should vectorize   on machines with   hardware  support  for */
/*       vector gather/scatter operations.  Your compiler may require */
/*       a  compiler directive  to  convince   it that there  are  no */
/*       implicit vector  dependencies.  Compiler directives  for the */
/*       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied */
/*       with the standard SLAP distribution. */

/* *Precision:           Single Precision */
/* *Cautions: */
/*     This   routine   assumes  that  the matrix A is stored in SLAP */
/*     Column format.  It does not check  for  this (for  speed)  and */
/*     evil, ugly, ornery and nasty things  will happen if the matrix */
/*     data  structure  is,  in fact, not SLAP Column.  Beware of the */
/*     wrong data structure!!! */

/* *See Also: */
/*       SSMTV */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SSMV */

/*         Zero out the result vector. */
/* ***FIRST EXECUTABLE STATEMENT  SSMV */
    /* Parameter adjustments */
    --y;
    --x;
    --a;
    --ja;
    --ia;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = 0.f;
/* L10: */
    }

/*         Multiply by A. */

/* VD$R NOCONCUR */
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	ibgn = ja[icol];
	iend = ja[icol + 1] - 1;
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
	i__2 = iend;
	for (i__ = ibgn; i__ <= i__2; ++i__) {
	    y[ia[i__]] += a[i__] * x[icol];
/* L20: */
	}
/* L30: */
    }

    if (*isym == 1) {

/*         The matrix is non-symmetric.  Need to get the other half in... */
/*         This loops assumes that the diagonal is the first entry in */
/*         each column. */

	i__1 = *n;
	for (irow = 1; irow <= i__1; ++irow) {
	    jbgn = ja[irow] + 1;
	    jend = ja[irow + 1] - 1;
	    if (jbgn > jend) {
		goto L50;
	    }
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		y[irow] += a[j] * x[ia[j]];
/* L40: */
	    }
L50:
	    ;
	}
    }
    return 0;
/* ------------- LAST LINE OF SSMV FOLLOWS ---------------------------- */
} /* ssmv_ */

/* DECK SSMTV */
/* Subroutine */ int ssmtv_(integer *n, real *x, real *y, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, ibgn, iend, jbgn, jend, icol, irow;

/* ***BEGIN PROLOGUE  SSMTV */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSMTV-S), */
/*             Matrix transpose Vector Multiply, Sparse */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP Column Format Sparse Matrix (transpose) Vector Prdt. */
/*            Routine to calculate the sparse matrix vector product: */
/*            Y = A'*X, where ' denotes transpose. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER  N, NELT, IA(NELT), JA(N+1), ISYM */
/*     REAL     X(N), Y(N), A(NELT) */

/*     CALL SSMTV(N, X, Y, NELT, IA, JA, A, ISYM ) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* X      :IN       Real X(N). */
/*         The vector that should be multiplied by the transpose of */
/*         the matrix. */
/* Y      :OUT      Real Y(N). */
/*         The product of the transpose of the matrix and the vector. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(N+1). */
/* A      :IN       Integer A(NELT). */
/*         These arrays should hold the matrix A in the SLAP Column */
/*         format.  See "Description", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */

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

/*       With  the SLAP  format  the "inner  loops" of  this  routine */
/*       should vectorize   on machines with   hardware  support  for */
/*       vector gather/scatter operations.  Your compiler may require */
/*       a  compiler directive  to  convince   it that there  are  no */
/*       implicit vector  dependencies.  Compiler directives  for the */
/*       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied */
/*       with the standard SLAP distribution. */

/* *Precision:           Single Precision */
/* *Cautions: */
/*     This   routine   assumes  that  the matrix A is stored in SLAP */
/*     Column format.  It does not check  for  this (for  speed)  and */
/*     evil, ugly, ornery and nasty things  will happen if the matrix */
/*     data  structure  is,  in fact, not SLAP Column.  Beware of the */
/*     wrong data structure!!! */

/* *See Also: */
/*       SSMV */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SSMTV */

/*         Zero out the result vector. */
/* ***FIRST EXECUTABLE STATEMENT  SSMTV */
    /* Parameter adjustments */
    --y;
    --x;
    --a;
    --ja;
    --ia;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	y[i__] = 0.f;
/* L10: */
    }

/*         Multiply by A-Transpose. */
/*         A-Transpose is stored by rows... */
/* VD$R NOCONCUR */
    i__1 = *n;
    for (irow = 1; irow <= i__1; ++irow) {
	ibgn = ja[irow];
	iend = ja[irow + 1] - 1;
/* VD$ ASSOC */
	i__2 = iend;
	for (i__ = ibgn; i__ <= i__2; ++i__) {
	    y[irow] += a[i__] * x[ia[i__]];
/* L20: */
	}
/* L30: */
    }

    if (*isym == 1) {

/*         The matrix is non-symmetric.  Need to get the other half in... */
/*         This loops assumes that the diagonal is the first entry in */
/*         each column. */

	i__1 = *n;
	for (icol = 1; icol <= i__1; ++icol) {
	    jbgn = ja[icol] + 1;
	    jend = ja[icol + 1] - 1;
	    if (jbgn > jend) {
		goto L50;
	    }
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		y[ia[j]] += a[j] * x[icol];
/* L40: */
	    }
L50:
	    ;
	}
    }
    return 0;
/* ------------- LAST LINE OF SSMTV FOLLOWS ---------------------------- */
} /* ssmtv_ */

/* DECK SSDI */
/* Subroutine */ int ssdi_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, real *rwork, 
	integer *iwork)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, locd;

/* ***BEGIN PROLOGUE  SSDI */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213  (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSDI-S), */
/*             Linear system solve, Sparse, Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Diagonal Matrix Vector Multiply. */
/*            Routine to calculate the product  X = DIAG*B, */
/*            where DIAG is a diagonal matrix. */
/* ***DESCRIPTION */
/* *Usage: */
/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Vector to multiply the diagonal by. */
/* X      :OUT      Real X(N). */
/*         Result of DIAG*B. */
/* NELT   :DUMMY    Integer. */
/*         Retained for compatibility with SLAP MSOLVE calling sequence. */
/* IA     :DUMMY    Integer IA(NELT). */
/*         Retained for compatibility with SLAP MSOLVE calling sequence. */
/* JA     :DUMMY    Integer JA(N+1). */
/*         Retained for compatibility with SLAP MSOLVE calling sequence. */
/*  A     :DUMMY    Real A(NELT). */
/*         Retained for compatibility with SLAP MSOLVE calling sequence. */
/* ISYM   :DUMMY    Integer. */
/*         Retained for compatibility with SLAP MSOLVE calling sequence. */
/* RWORK  :IN       Real RWORK(USER DEFINABLE). */
/*         Work array holding the diagonal of some matrix to scale */
/*         B by.  This array must be set by the user or by a call */
/*         to the slap routine SSDS or SSD2S.  The length of RWORK */
/*         must be > IWORK(4)+N. */
/* IWORK  :IN       Integer IWORK(10). */
/*         IWORK(4) holds the offset into RWORK for the diagonal matrix */
/*         to scale B by.  This is usually set up by the SLAP pre- */
/*         conditioner setup routines SSDS or SSD2S. */

/* *Description: */
/*         This routine is supplied with the SLAP package to perform */
/*         the  MSOLVE  operation for iterative drivers that require */
/*         diagonal  Scaling  (e.g., SSDCG, SSDBCG).   It  conforms */
/*         to the SLAP MSOLVE CALLING CONVENTION  and hence does not */
/*         require an interface routine as do some of the other pre- */
/*         conditioners supplied with SLAP. */

/* *Precision:           Single Precision */
/* *See Also: */
/*       SSDS, SSD2S */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SSDI */

/*         Determine where the inverse of the diagonal */
/*         is in the work array and then scale by it. */
/* ***FIRST EXECUTABLE STATEMENT  SSDI */
    /* Parameter adjustments */
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    locd = iwork[4] - 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = rwork[locd + i__] * b[i__];
/* L10: */
    }
    return 0;
/* ------------- LAST LINE OF SSDI FOLLOWS ---------------------------- */
} /* ssdi_ */

/* DECK SSLI */
/* Subroutine */ int ssli_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, real *rwork, 
	integer *iwork)
{
    static integer nel;
    extern /* Subroutine */ int ssli2_(integer *, real *, real *, integer *, 
	    integer *, integer *, real *);
    static integer locel, lociel, locjel;

/* ***BEGIN PROLOGUE  SSLI */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSLI-S), */
/*             Linear system solve, Sparse, Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP MSOLVE for Lower Triangle Matrix. */
/*            This routine acts as an interface between the SLAP generic */
/*            MSOLVE calling convention and the routine that actually */
/*                      -1 */
/*            computes L  B = X. */

/* *Description */
/*       See the Description of SLLI2 for the gory details. */
/* ***ROUTINES CALLED  SLLI2 */
/* ***END PROLOGUE  SSLI */
/* ***FIRST EXECUTABLE STATEMENT  SSLI */

    /* Parameter adjustments */
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    nel = iwork[1];
    lociel = iwork[2];
    locjel = iwork[3];
    locel = iwork[4];
    ssli2_(n, &b[1], &x[1], &nel, &iwork[lociel], &iwork[locjel], &rwork[
	    locel]);

    return 0;
/* ------------- LAST LINE OF SSLI FOLLOWS ---------------------------- */
} /* ssli_ */

/* DECK SSLI2 */
/* Subroutine */ int ssli2_(integer *n, real *b, real *x, integer *nel, 
	integer *iel, integer *jel, real *el)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, jbgn, jend, icol;

/* ***BEGIN PROLOGUE  SSLI2 */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSLI2-S), */
/*             Linear system solve, Sparse, Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP for Lower Triangle Matrix Backsolve. */
/*            Routine to solve a system of the form  Lx = b , where */
/*            L is a lower triangular matrix. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N,  NEL, IEL(N+1), JEL(NEL) */
/*     REAL    B(N), X(N), EL(NEL) */

/*     CALL SSLI2( N, B, X, NEL, IEL, JEL, EL ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right hand side vector. */
/* X      :OUT      Real X(N). */
/*         Solution to Lx = b. */
/* NEL    :IN       Integer. */
/*         Number of non-zeros in the EL array. */
/* IEL    :IN       Integer IEL(N+1). */
/* JEL    :IN       Integer JEL(NEL). */
/* EL     :IN       Real     EL(NEL). */
/*         IEL, JEL, EL contain the unit lower triangular factor   of */
/*         the incomplete decomposition   of the A  matrix  stored in */
/*         SLAP Row format.  The diagonal of  ones *IS* stored.  This */
/*         structure can be set up by the  SS2LT  routine.  See "LONG */
/*         DESCRIPTION", below for more details about  the  SLAP  Row */
/*         format. */

/* *Description: */
/*       This routine is supplied with the SLAP package  as a routine */
/*       to  perform the  MSOLVE operation in  the SIR for the driver */
/*       routine SSGS.  It must be called via the SLAP MSOLVE calling */
/*       sequence convention interface routine SSLI. */
/*         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE **** */
/*               **** SLAP MSOLVE CALLING CONVENTION **** */

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

/*       With  the SLAP  Row format  the "inner loop" of this routine */
/*       should vectorize   on machines with   hardware  support  for */
/*       vector gather/scatter operations.  Your compiler may require */
/*       a  compiler directive  to  convince   it that there  are  no */
/*       implicit vector  dependencies.  Compiler directives  for the */
/*       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied */
/*       with the standard SLAP distribution. */

/* *Precision: Single Precision */
/* *See Also: */
/*         SSLI */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SSLI2 */

/*         Initialize the solution by copying the right hands side */
/*         into it. */
/* ***FIRST EXECUTABLE STATEMENT  SSLI2 */
    /* Parameter adjustments */
    --x;
    --b;
    --el;
    --jel;
    --iel;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = b[i__];
/* L10: */
    }

/* VD$ NOCONCUR */
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	x[icol] /= el[jel[icol]];
	jbgn = jel[icol] + 1;
	jend = jel[icol + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NOCONCUR */
/* VD$ NODEPCHK */
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		x[iel[j]] -= el[j] * x[icol];
/* L20: */
	    }
	}
/* L30: */
    }

    return 0;
/* ------------- LAST LINE OF SSLI2 FOLLOWS ---------------------------- */
} /* ssli2_ */

/* DECK SSLLTI */
/* Subroutine */ int ssllti_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, real *rwork, 
	integer *iwork)
{
    static integer nel, locel;
    extern /* Subroutine */ int sllti2_(integer *, real *, real *, integer *, 
	    integer *, integer *, real *, real *);
    static integer lociel, locjel, locdin;

/* ***BEGIN PROLOGUE  SSLLTI */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSLLTI-S), */
/*             Linear system solve, Sparse, Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP MSOLVE for LDL' (IC) Factorization. */
/*            This routine acts as an interface between the SLAP generic */
/*            MSOLVE calling convention and the routine that actually */
/*                           -1 */
/*            computes (LDL')  B = X. */
/* ***DESCRIPTION */
/*       See the DESCRIPTION of SLLTI2 for the gory details. */
/* ***ROUTINES CALLED SLLTI2 */

/* ***END PROLOGUE  SSLLTI */

/* ***FIRST EXECUTABLE STATEMENT  SSLLTI */
    /* Parameter adjustments */
    --b;
    --x;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    nel = iwork[1];
    lociel = iwork[3];
    locjel = iwork[2];
    locel = iwork[4];
    locdin = iwork[5];
    sllti2_(n, &b[1], &x[1], &nel, &iwork[lociel], &iwork[locjel], &rwork[
	    locel], &rwork[locdin]);

    return 0;
/* ------------- LAST LINE OF SSLLTI FOLLOWS ---------------------------- */
} /* ssllti_ */

/* DECK SLLTI2 */
/* Subroutine */ int sllti2_(integer *n, real *b, real *x, integer *nel, 
	integer *iel, integer *jel, real *el, real *dinv)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, ibgn, iend, irow;

/* ***BEGIN PROLOGUE  SLLTI2 */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SLLTI2-S), */
/*             Symmetric Linear system solve, Sparse, */
/*             Iterative Precondition, Incomplete Factorization */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP back solve routine for LDL' Factorization. */
/*            Routine to solve a system of the  form  L*D*L' X  =  B, */
/*            where L is a unit lower triangular  matrix  and  D is a */
/*            diagonal matrix and ' means transpose. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N,  NEL, IEL(N+1), JEL(NEL) */
/*     REAL    B(N), X(N), EL(NEL), DINV(N) */

/*     CALL SLLTI2( N, B, X, NEL, IEL, JEL, EL, DINV ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right hand side vector. */
/* X      :OUT      Real X(N). */
/*         Solution to L*D*L' x = b. */
/* NEL    :IN       Integer. */
/*         Number of non-zeros in the EL array. */
/* IEL    :IN       Integer IEL(N+1). */
/* JEL    :IN       Integer JEL(NEL). */
/* EL     :IN       Real     EL(NEL). */
/*         IEL, JEL, EL contain the unit lower triangular factor   of */
/*         the incomplete decomposition   of the A  matrix  stored in */
/*         SLAP Row format.   The diagonal of ones *IS* stored.  This */
/*         structure can be set  up  by  the SS2LT routine. See */
/*         "Description", below for more details about the   SLAP Row */
/*         format. */
/* DINV   :IN       Real DINV(N). */
/*         Inverse of the diagonal matrix D. */

/* *Description: */
/*       This routine is supplied with  the SLAP package as a routine */
/*       to perform the MSOLVE operation in the SCG iteration routine */
/*       for  the driver  routine SSICCG.   It must be called via the */
/*       SLAP  MSOLVE calling sequence  convention  interface routine */
/*       SSLLI. */
/*         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE **** */
/*               **** SLAP MSOLVE CALLING CONVENTION **** */

/*       IEL, JEL, EL should contain the unit lower triangular factor */
/*       of  the incomplete decomposition of  the A matrix  stored in */
/*       SLAP Row format.   This IC factorization  can be computed by */
/*       the  SSICS routine.  The  diagonal  (which is all one's) is */
/*       stored. */

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

/*       With  the SLAP  Row format  the "inner loop" of this routine */
/*       should vectorize   on machines with   hardware  support  for */
/*       vector gather/scatter operations.  Your compiler may require */
/*       a  compiler directive  to  convince   it that there  are  no */
/*       implicit vector  dependencies.  Compiler directives  for the */
/*       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied */
/*       with the standard SLAP distribution. */

/* *Precision:           Single Precision */
/* *See Also: */
/*       SSICCG, SSICS */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SLLTI2 */

/*         solve  l*y = b,  storing result in x. */
/* ***FIRST EXECUTABLE STATEMENT  SLLTI2 */
    /* Parameter adjustments */
    --dinv;
    --x;
    --b;
    --el;
    --iel;
    --jel;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = b[i__];
/* L10: */
    }
    i__1 = *n;
    for (irow = 1; irow <= i__1; ++irow) {
	ibgn = iel[irow] + 1;
	iend = iel[irow + 1] - 1;
	if (ibgn <= iend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NOCONCUR */
/* VD$ NODEPCHK */
	    i__2 = iend;
	    for (i__ = ibgn; i__ <= i__2; ++i__) {
		x[irow] -= el[i__] * x[jel[i__]];
/* L20: */
	    }
	}
/* L30: */
    }

/*         Solve  D*Z = Y,  storing result in X. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] *= dinv[i__];
/* L40: */
    }

/*         Solve  L-trans*X = Z. */

    for (irow = *n; irow >= 2; --irow) {
	ibgn = iel[irow] + 1;
	iend = iel[irow + 1] - 1;
	if (ibgn <= iend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NOCONCUR */
/* VD$ NODEPCHK */
	    i__1 = iend;
	    for (i__ = ibgn; i__ <= i__1; ++i__) {
		x[jel[i__]] -= el[i__] * x[irow];
/* L50: */
	    }
	}
/* L60: */
    }

    return 0;
/* ------------- LAST LINE OF SLTI2 FOLLOWS ---------------------------- */
} /* sllti2_ */

/* DECK SSLUI */
/* Subroutine */ int sslui_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, real *rwork, 
	integer *iwork)
{
    static integer locl, locu, locil, locjl, lociu, locju;
    extern /* Subroutine */ int sslui2_(integer *, real *, real *, integer *, 
	    integer *, real *, real *, integer *, integer *, real *);
    static integer locdin;

/* ***BEGIN PROLOGUE  SSLUI */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSLUI-S), */
/*             Non-Symmetric Linear system solve, Sparse, */
/*             Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP MSOLVE for LDU Factorization. */
/*            This routine  acts as an  interface between  the   SLAP */
/*            generic MSLOVE calling convention and the routine  that */
/*            actually computes:     -1 */
/*                              (LDU)  B = X. */
/* ***DESCRIPTION */
/*       See the "DESCRIPTION" of SSLUI2 for the gory details. */
/* ***ROUTINES CALLED  SSLUI2 */
/* ***END PROLOGUE  SSLUI */

/*         Pull out the locations of the arrays holding the ILU */
/*         factorization. */
/* ***FIRST EXECUTABLE STATEMENT  SSLUI */
    /* Parameter adjustments */
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    locil = iwork[1];
    locjl = iwork[2];
    lociu = iwork[3];
    locju = iwork[4];
    locl = iwork[5];
    locdin = iwork[6];
    locu = iwork[7];

/*         Solve the system LUx = b */
    sslui2_(n, &b[1], &x[1], &iwork[locil], &iwork[locjl], &rwork[locl], &
	    rwork[locdin], &iwork[lociu], &iwork[locju], &rwork[locu]);

    return 0;
/* ------------- LAST LINE OF SSLUI FOLLOWS ---------------------------- */
} /* sslui_ */

/* DECK SSLUI2 */
/* Subroutine */ int sslui2_(integer *n, real *b, real *x, integer *il, 
	integer *jl, real *l, real *dinv, integer *iu, integer *ju, real *u)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, jbgn, jend, icol, irow;

/* ***BEGIN PROLOGUE  SSLUI2 */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSLUI2-S), */
/*             Non-Symmetric Linear system solve, Sparse, */
/*             Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP Back solve for LDU Factorization. */
/*            Routine  to  solve a system of the form  L*D*U X  =  B, */
/*            where L is a unit  lower  triangular  matrix,  D  is  a */
/*            diagonal matrix, and U is a unit upper triangular matrix. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, IL(N+1), JL(NL), IU(NU), JU(N+1) */
/*     REAL    B(N), X(N), L(NL), DINV(N), U(NU) */

/*     CALL SSLUI2( N, B, X, IL, JL, L, DINV, IU, JU, U ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right hand side. */
/* X      :OUT      Real X(N). */
/*         Solution of L*D*U x = b. */
/* NEL    :IN       Integer. */
/*         Number of non-zeros in the EL array. */
/* IL     :IN       Integer IL(N+1). */
/* JL     :IN       Integer JL(NL). */
/*  L     :IN       Real     L(NL). */
/*         IL, JL, L contain the unit  lower triangular factor of the */
/*         incomplete decomposition of some matrix stored in SLAP Row */
/*         format.  The diagonal of ones *IS* stored.  This structure */
/*         can   be   set up  by   the  SSILUS routine.   See */
/*         "DESCRIPTION", below  for more   details about   the  SLAP */
/*         format. */
/* DINV   :IN       Real DINV(N). */
/*         Inverse of the diagonal matrix D. */
/* NU     :IN       Integer. */
/*         Number of non-zeros in the U array. */
/* IU     :IN       Integer IU(N+1). */
/* JU     :IN       Integer JU(NU). */
/* U      :IN       Real     U(NU). */
/*         IU, JU, U contain the unit upper triangular factor  of the */
/*         incomplete decomposition  of  some  matrix stored in  SLAP */
/*         Column format.   The diagonal of ones  *IS* stored.   This */
/*         structure can be set up  by the SSILUS routine.  See */
/*         "DESCRIPTION", below   for  more   details about  the SLAP */
/*         format. */

/* *Description: */
/*       This routine is supplied with  the SLAP package as a routine */
/*       to  perform  the  MSOLVE operation  in   the  SIR and   SBCG */
/*       iteration routines for  the  drivers SSILUR and SSLUBC.   It */
/*       must  be called  via   the  SLAP  MSOLVE  calling   sequence */
/*       convention interface routine SSLUI. */
/*         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE **** */
/*               **** SLAP MSOLVE CALLING CONVENTION **** */

/*       IL, JL, L should contain the unit lower triangular factor of */
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

/*       With  the SLAP  format  the "inner  loops" of  this  routine */
/*       should vectorize   on machines with   hardware  support  for */
/*       vector gather/scatter operations.  Your compiler may require */
/*       a  compiler directive  to  convince   it that there  are  no */
/*       implicit vector  dependencies.  Compiler directives  for the */
/*       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied */
/*       with the standard SLAP distribution. */

/* *Precision:           Single Precision */
/* *See Also: */
/*       SSILUS */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SSLUI2 */

/*         Solve  L*Y = B,  storing result in X, L stored by rows. */
/* ***FIRST EXECUTABLE STATEMENT  SSLUI2 */
    /* Parameter adjustments */
    --dinv;
    --x;
    --b;
    --il;
    --jl;
    --l;
    --iu;
    --ju;
    --u;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = b[i__];
/* L10: */
    }
    i__1 = *n;
    for (irow = 2; irow <= i__1; ++irow) {
	jbgn = il[irow];
	jend = il[irow + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ ASSOC */
/* VD$ NODEPCHK */
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		x[irow] -= l[j] * x[jl[j]];
/* L20: */
	    }
	}
/* L30: */
    }

/*         Solve  D*Z = Y,  storing result in X. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] *= dinv[i__];
/* L40: */
    }

/*         Solve  U*X = Z, U stored by columns. */
    for (icol = *n; icol >= 2; --icol) {
	jbgn = ju[icol];
	jend = ju[icol + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
	    i__1 = jend;
	    for (j = jbgn; j <= i__1; ++j) {
		x[iu[j]] -= u[j] * x[icol];
/* L50: */
	    }
	}
/* L60: */
    }

    return 0;
/* ------------- LAST LINE OF SSLUI2 FOLLOWS ---------------------------- */
} /* sslui2_ */

/* DECK SSLUTI */
/* Subroutine */ int ssluti_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, real *rwork, 
	integer *iwork)
{
    static integer locl, locu, locil, locjl, lociu, locju;
    extern /* Subroutine */ int sslui4_(integer *, real *, real *, integer *, 
	    integer *, real *, real *, integer *, integer *, real *);
    static integer locdin;

/* ***BEGIN PROLOGUE  SSLUTI */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSLUTI-S), */
/*             Linear system solve, Sparse, Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP MTSOLV for LDU Factorization. */
/*            This routine acts as  an  interface  between  the  SLAP */
/*            generic MTSOLV calling convention and  the routine that */
/*            actually computes:       -T */
/*                                (LDU)  B = X. */
/* ***DESCRIPTION */
/*       See the "DESCRIPTION" of SSLUI4 for the gory details. */
/* ***ROUTINES CALLED  SSLUI4 */
/* ***END PROLOGUE  SSLUTI */

/*         Pull out the pointers to the L, D and U matricies and call */
/*         the workhorse routine. */
/* ***FIRST EXECUTABLE STATEMENT  SSLUTI */
    /* Parameter adjustments */
    --a;
    --x;
    --b;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    locil = iwork[1];
    locjl = iwork[2];
    lociu = iwork[3];
    locju = iwork[4];
    locl = iwork[5];
    locdin = iwork[6];
    locu = iwork[7];

    sslui4_(n, &b[1], &x[1], &iwork[locil], &iwork[locjl], &rwork[locl], &
	    rwork[locdin], &iwork[lociu], &iwork[locju], &rwork[locu]);

    return 0;
/* ------------- LAST LINE OF SSLUTI FOLLOWS ---------------------------- */
} /* ssluti_ */

/* DECK SSLUI4 */
/* Subroutine */ int sslui4_(integer *n, real *b, real *x, integer *il, 
	integer *jl, real *l, real *dinv, integer *iu, integer *ju, real *u)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, jbgn, jend, icol, irow;

/* ***BEGIN PROLOGUE  SSLUI4 */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSLUI4-S), */
/*             Non-Symmetric Linear system solve, Sparse, */
/*             Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP back solve for LDU Factorization. */
/*            Routine to solve a system of the form  (L*D*U)' X =  B, */
/*            where L is a unit  lower  triangular  matrix,  D  is  a */
/*            diagonal matrix, and  U  is  a  unit  upper  triangular */
/*            matrix and ' denotes transpose. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NL, IL(N+1), JL(NL), NU, IU(N+1), JU(NU) */
/*     REAL    B(N), X(N), L(NEL), DINV(N), U(NU) */

/*     CALL SSLUI4( N, B, X, IL, JL, L, DINV, IU, JU, U ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right hand side. */
/* X      :OUT      Real X(N). */
/*         Solution of (L*D*U)trans x = b. */
/* IL     :IN       Integer IL(N+1). */
/* JL     :IN       Integer JL(NL). */
/*  L     :IN       Real     L(NL). */
/*         IL, JL, L contain the unit lower triangular  factor of the */
/*         incomplete decomposition of some matrix stored in SLAP Row */
/*         format.  The diagonal of ones *IS* stored.  This structure */
/*         can    be set  up  by   the  SSILUS routine.   See */
/*         "DESCRIPTION",  below for  more  details about  the   SLAP */
/*         format. */
/* DINV   :IN       Real DINV(N). */
/*         Inverse of the diagonal matrix D. */
/* IU     :IN       Integer IU(N+1). */
/* JU     :IN       Integer JU(NU). */
/* U      :IN       Real     U(NU). */
/*         IU, JU, U contain the  unit upper triangular factor of the */
/*         incomplete  decomposition of some  matrix stored  in  SLAP */
/*         Column  format.   The diagonal of  ones *IS* stored.  This */
/*         structure can be set up by the  SSILUS routine.  See */
/*         "DESCRIPTION",  below for  more  details  about  the  SLAP */
/*         format. */

/* *Description: */
/*       This routine is supplied with the SLAP package as  a routine */
/*       to  perform  the  MTSOLV  operation  in  the SBCG  iteration */
/*       routine for the  driver SSLUBC.   It must  be called via the */
/*       SLAP  MTSOLV calling  sequence convention interface  routine */
/*       SSLUTI. */
/*         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE **** */
/*               **** SLAP MSOLVE CALLING CONVENTION **** */

/*       IL, JL, L should contain the unit lower triangular factor of */
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

/*       With  the SLAP  format  the "inner  loops" of  this  routine */
/*       should vectorize   on machines with   hardware  support  for */
/*       vector gather/scatter operations.  Your compiler may require */
/*       a  compiler directive  to  convince   it that there  are  no */
/*       implicit vector  dependencies.  Compiler directives  for the */
/*       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied */
/*       with the standard SLAP distribution. */

/* *Precision:           Single Precision */
/* *See Also: */
/*       SSILUS */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SSLUI4 */

/* ***FIRST EXECUTABLE STATEMENT  SSLUI4 */
    /* Parameter adjustments */
    --dinv;
    --x;
    --b;
    --il;
    --jl;
    --l;
    --iu;
    --ju;
    --u;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = b[i__];
/* L10: */
    }

/*         Solve  U'*Y = X,  storing result in X, U stored by columns. */
    i__1 = *n;
    for (irow = 2; irow <= i__1; ++irow) {
	jbgn = ju[irow];
	jend = ju[irow + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ ASSOC */
/* VD$ NODEPCHK */
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		x[irow] -= u[j] * x[iu[j]];
/* L70: */
	    }
	}
/* L80: */
    }

/*         Solve  D*Z = Y,  storing result in X. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] *= dinv[i__];
/* L90: */
    }

/*         Solve  L'*X = Z, L stored by rows. */
    for (icol = *n; icol >= 2; --icol) {
	jbgn = il[icol];
	jend = il[icol + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
	    i__1 = jend;
	    for (j = jbgn; j <= i__1; ++j) {
		x[jl[j]] -= l[j] * x[icol];
/* L100: */
	    }
	}
/* L110: */
    }
    return 0;
/* ------------- LAST LINE OF SSLUI4 FOLLOWS ---------------------------- */
} /* sslui4_ */

/* DECK SSMMTI */
/* Subroutine */ int ssmmti_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, real *rwork, 
	integer *iwork)
{
    static integer locl, locu, locil, locjl, lociu, locju;
    extern /* Subroutine */ int ssmmi2_(integer *, real *, real *, integer *, 
	    integer *, real *, real *, integer *, integer *, real *);
    static integer locdin;

/* ***BEGIN PROLOGUE  SSMMTI */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSMMTI-S), */
/*             Linear system solve, Sparse, Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP MSOLVE for LDU Factorization of Normal Equations. */
/*            This routine acts as  an  interface  between  the  SLAP */
/*            generic MMTSLV calling convention and the routine  that */
/*            actually computes:            -1 */
/*                            [(LDU)*(LDU)']  B = X. */
/* ***DESCRIPTION */
/*       See the "DESCRIPTION" of SSMMI2 for the gory details. */
/* ***ROUTINES CALLED  SSMMI2 */
/* ***END PROLOGUE  SSMMTI */

/*         Pull out the locations of the arrays holding the ILU */
/*         factorization. */
/* ***FIRST EXECUTABLE STATEMENT  SSMMTI */
    /* Parameter adjustments */
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    locil = iwork[1];
    locjl = iwork[2];
    lociu = iwork[3];
    locju = iwork[4];
    locl = iwork[5];
    locdin = iwork[6];
    locu = iwork[7];

    ssmmi2_(n, &b[1], &x[1], &iwork[locil], &iwork[locjl], &rwork[locl], &
	    rwork[locdin], &iwork[lociu], &iwork[locju], &rwork[locu]);

    return 0;
/* ------------- LAST LINE OF SSMMTI FOLLOWS ---------------------------- */
} /* ssmmti_ */

/* DECK SSMMI2 */
/* Subroutine */ int ssmmi2_(integer *n, real *b, real *x, integer *il, 
	integer *jl, real *l, real *dinv, integer *iu, integer *ju, real *u)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, jbgn, jend, icol, irow;

/* ***BEGIN PROLOGUE  SSMMI2 */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSMMI2-S), */
/*             Linear system, Sparse, Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP Back solve for LDU Factorization of Normal Equations. */
/*            To solve a system of the form (L*D*U)*(L*D*U)' X  =  B, */
/*            where  L  is a unit lower triangular matrix,  D   is  a */
/*            diagonal matrix, and  U  is  a  unit  upper  triangular */
/*            matrix and ' denotes transpose. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, IL(N+1), JL(NL), IU(N+1), JU(NU) */
/*     REAL    B(N), X(N), L(NL), DINV(N), U(NU) */

/*     CALL SSMMI2( N, B, X, IL, JE, L, DINV, IU, JU, U ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right hand side. */
/* X      :OUT      Real X(N). */
/*         Solution of (L*D*U)(L*D*U)trans x = b. */
/* IL     :IN       Integer IL(N+1). */
/* JL     :IN       Integer JL(NL). */
/*  L     :IN       Real     L(NL). */
/*         IL, JL, L contain the unit lower  triangular factor of the */
/*         incomplete decomposition of some matrix stored in SLAP Row */
/*         format.  The diagonal of ones *IS* stored.  This structure */
/*         can  be  set up by   the  SSILUS   routine.    See */
/*         "DESCRIPTION", below for  more   details   about  the SLAP */
/*         format. */
/* DINV   :IN       Real DINV(N). */
/*         Inverse of the diagonal matrix D. */
/* IU     :IN       Integer IU(N+1). */
/* JU     :IN       Integer JU(NU). */
/* U      :IN       Real     U(NU). */
/*         IU, JU, U contain the unit upper  triangular factor of the */
/*         incomplete decomposition  of   some matrix stored in  SLAP */
/*         Column  format.  The diagonal  of  ones *IS* stored.  This */
/*         structure can be set up  by the SSILUS routine.  See */
/*         "DESCRIPTION",  below  for  more  details  about  the SLAP */
/*         format. */

/* *Description: */
/*       This routine is supplied with the SLAP package as  a routine */
/*       to  perform  the  MSOLVE  operation  in  the SBCGN iteration */
/*       routine for the  driver SSLUCN.   It must  be called via the */
/*       SLAP  MSOLVE calling  sequence convention interface  routine */
/*       SSMMTI. */
/*         **** THIS ROUTINE ITSELF DOES NOT CONFORM TO THE **** */
/*               **** SLAP MSOLVE CALLING CONVENTION **** */

/*       IL, JL, L should contain the unit lower triangular factor of */
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

/*       With  the SLAP  format  the "inner  loops" of  this  routine */
/*       should vectorize   on machines with   hardware  support  for */
/*       vector gather/scatter operations.  Your compiler may require */
/*       a  compiler directive  to  convince   it that there  are  no */
/*       implicit vector  dependencies.  Compiler directives  for the */
/*       Alliant FX/Fortran and CRI CFT/CFT77 compilers  are supplied */
/*       with the standard SLAP distribution. */

/* *Precision:           Single Precision */
/* *See Also: */
/*       SSILUS */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SSMMI2 */

/*         Solve  L*Y = B,  storing result in X, L stored by rows. */
/* ***FIRST EXECUTABLE STATEMENT  SSMMI2 */
    /* Parameter adjustments */
    --u;
    --dinv;
    --x;
    --b;
    --il;
    --jl;
    --l;
    --iu;
    --ju;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = b[i__];
/* L10: */
    }
    i__1 = *n;
    for (irow = 2; irow <= i__1; ++irow) {
	jbgn = il[irow];
	jend = il[irow + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ ASSOC */
/* VD$ NODEPCHK */
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		x[irow] -= l[j] * x[jl[j]];
/* L20: */
	    }
	}
/* L30: */
    }

/*         Solve  D*Z = Y,  storing result in X. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] *= dinv[i__];
/* L40: */
    }

/*         Solve  U*X = Z, U stored by columns. */
    for (icol = *n; icol >= 2; --icol) {
	jbgn = ju[icol];
	jend = ju[icol + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
	    i__1 = jend;
	    for (j = jbgn; j <= i__1; ++j) {
		x[iu[j]] -= u[j] * x[icol];
/* L50: */
	    }
	}
/* L60: */
    }

/*         Solve  U'*Y = X,  storing result in X, U stored by columns. */
    i__1 = *n;
    for (irow = 2; irow <= i__1; ++irow) {
	jbgn = ju[irow];
	jend = ju[irow + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ ASSOC */
/* VD$ NODEPCHK */
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		x[irow] -= u[j] * x[iu[j]];
/* L70: */
	    }
	}
/* L80: */
    }

/*         Solve  D*Z = Y,  storing result in X. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] *= dinv[i__];
/* L90: */
    }

/*         Solve  L'*X = Z, L stored by rows. */
    for (icol = *n; icol >= 2; --icol) {
	jbgn = il[icol];
	jend = il[icol + 1] - 1;
	if (jbgn <= jend) {
/* LLL. OPTION ASSERT (NOHAZARD) */
/* DIR$ IVDEP */
/* VD$ NODEPCHK */
	    i__1 = jend;
	    for (j = jbgn; j <= i__1; ++j) {
		x[jl[j]] -= l[j] * x[icol];
/* L100: */
	    }
	}
/* L110: */
    }

    return 0;
/* ------------- LAST LINE OF SSMMI2 FOLLOWS ---------------------------- */
} /* ssmmi2_ */

