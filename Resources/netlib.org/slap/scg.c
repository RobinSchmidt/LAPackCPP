/* scg.f -- translated by f2c (version 20100827).
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

static integer c__3 = 3;
static integer c__1 = 1;
static integer c__113 = 113;
static integer c__0 = 0;
static real c_b28 = 0.f;

/* DECK SCG */
/* Subroutine */ int scg_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, S_fp matvec, S_fp 
	msolve, integer *itol, real *tol, integer *itmax, integer *iter, real 
	*err, integer *ierr, integer *iunit, real *r__, real *z__, real *p, 
	real *dz, real *rwork, integer *iwork)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, k;
    static real ak, bk, bnrm;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real akden, bkden;
    extern integer isscg_(integer *, real *, real *, integer *, integer *, 
	    integer *, real *, integer *, S_fp, integer *, real *, integer *, 
	    integer *, real *, integer *, integer *, real *, real *, real *, 
	    real *, real *, integer *, real *, real *, real *, real *);
    static real bknum;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), saxpy_(integer *, real *, real *, integer *, real *, 
	    integer *);
    extern doublereal r1mach_(integer *);
    static real tolmin, solnrm;

/* ***BEGIN PROLOGUE  SCG */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SCG-S), */
/*             Symmetric Linear system, Sparse, Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Preconditioned Conjugate Gradient iterative Ax=b solver. */
/*            Routine to  solve a  symmetric positive definite linear */
/*            system    Ax = b    using the Preconditioned  Conjugate */
/*            Gradient method. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX */
/*     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINABLE) */
/*     REAL     B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N), P(N), DZ(N) */
/*     REAL     RWORK(USER DEFINABLE) */
/*     EXTERNAL MATVEC, MSOLVE */

/*     CALL SCG(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSLOVE, */
/*    $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, */
/*    $     RWORK, IWORK ) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right-hand side vector. */
/* X      :INOUT    Real X(N). */
/*         On input X is your initial guess for solution vector. */
/*         On output X is the final approximate solution. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Integer A(NELT). */
/*         These arrays contain the matrix data structure for A. */
/*         It could take any form.  See ``Description'', below */
/*         for more late breaking details... */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* MATVEC :EXT      External. */
/*         Name of a routine which performs the matrix vector multiply */
/*         Y = A*X given A and X.  The name of the MATVEC routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MATVEC is: */

/*             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM ) */

/*         Where N is the number of unknowns, Y is the product A*X */
/*         upon return X is an input vector, NELT is the number of */
/*         non-zeros in the SLAP IA, JA, A storage for the matrix A. */
/*         ISYM is a flag which, if non-zero, denotest that A is */
/*         symmetric and only the lower or upper triangle is stored. */
/* MSOLVE :EXT      External. */
/*         Name of a routine which solves a linear system MZ = R for */
/*         Z given R with the preconditioning matrix M (M is supplied via */
/*         RWORK and IWORK arrays).  The name of the MSOLVE routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MSLOVE is: */

/*             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK) */

/*         Where N is the number of unknowns, R is the right-hand side */
/*         vector, and Z is the solution upon return.  RWORK is a real */
/*         array that can be used to pass necessary preconditioning */
/*         information and/or workspace to MSOLVE.  IWORK is an integer */
/*         work array for the same purpose as RWORK. */
/* ITOL   :IN       Integer. */
/*         Flag to indicate type of convergence criterion. */
/*         If ITOL=1, iteration stops when the 2-norm of the residual */
/*         divided by the 2-norm of the right-hand side is less than TOL. */
/*         If ITOL=2, iteration stops when the 2-norm of M-inv times the */
/*         residual divided by the 2-norm of M-inv times the right hand */
/*         side is less than TOL, where M-inv is the inverse of the */
/*         diagonal of A. */
/*         ITOL=11 is often useful for checking and comparing different */
/*         routines.  For this case, the user must supply the "exact" */
/*         solution or a very accurate approximation (one with an error */
/*         much less than TOL) through a common block, */
/*                     COMMON /SOLBLK/ SOLN(1) */
/*         if ITOL=11, iteration stops when the 2-norm of the difference */
/*         between the iterative approximation and the user-supplied */
/*         solution divided by the 2-norm of the user-supplied solution */
/*         is less than TOL.  Note that this requires the user to set up */
/*         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine. */
/*         The routine with this declaration should be loaded before the */
/*         stop test so that the correct length is used by the loader. */
/*         This procedure is not standard Fortran and may not work */
/*         correctly on your system (although it has worked on every */
/*         system the authors have tried).  If ITOL is not 11 then this */
/*         common block is indeed standard Fortran. */
/* TOL    :IN       Real. */
/*         Convergence criterion, as described above. */
/* ITMAX  :IN       Integer. */
/*         Maximum number of iterations. */
/* ITER   :OUT      Integer. */
/*         Number of iterations required to reach convergence, or */
/*         ITMAX+1 if convergence criterion could not be achieved in */
/*         ITMAX iterations. */
/* ERR    :OUT      Real. */
/*         Error estimate of error in final approximate solution, as */
/*         defined by ITOL. */
/* IERR   :OUT      Integer. */
/*         Return error flag. */
/*           IERR = 0 => All went well. */
/*           IERR = 1 => Insufficient storage allocated */
/*                       for WORK or IWORK. */
/*           IERR = 2 => Method failed to converge in */
/*                       ITMAX steps. */
/*           IERR = 3 => Error in user input.  Check input */
/*                       value of N, ITOL. */
/*           IERR = 4 => User error tolerance set too tight. */
/*                       Reset to 500.0*R1MACH(3).  Iteration proceeded. */
/*           IERR = 5 => Preconditioning matrix, M,  is not */
/*                       Positive Definite.  $(r,z) < 0.0$. */
/*           IERR = 6 => Matrix A is not Positive Definite. */
/*                       $(p,Ap) < 0.0$. */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* R      :WORK     Real R(N). */
/* Z      :WORK     Real Z(N). */
/* P      :WORK     Real P(N). */
/* DZ     :WORK     Real DZ(N). */
/* RWORK  :WORK     Real RWORK(USER DEFINABLE). */
/*         Real array that can be used by  MSOLVE. */
/* IWORK  :WORK     Integer IWORK(USER DEFINABLE). */
/*         Integer array that can be used by  MSOLVE. */

/* *Description */
/*       This routine does  not care  what matrix data   structure is */
/*       used for  A and M.  It simply   calls  the MATVEC and MSOLVE */
/*       routines, with  the arguments as  described above.  The user */
/*       could write any type of structure and the appropriate MATVEC */
/*       and MSOLVE routines.  It is assumed  that A is stored in the */
/*       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is */
/*       stored  in  IWORK  and  RWORK   in  some fashion.   The SLAP */
/*       routines SSDCG and SSICCG are examples of this procedure. */

/*       Two  examples  of  matrix  data structures  are the: 1) SLAP */
/*       Triad  format and 2) SLAP Column format. */

/*       =================== S L A P Triad format =================== */

/*       In  this   format only the  non-zeros are  stored.  They may */
/*       appear  in *ANY* order.   The user  supplies three arrays of */
/*       length NELT, where  NELT  is the number  of non-zeros in the */
/*       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero */
/*       the  user puts   the row  and  column index   of that matrix */
/*       element in the IA and JA arrays.  The  value of the non-zero */
/*       matrix  element is  placed in  the corresponding location of */
/*       the A  array.  This is  an extremely easy data  structure to */
/*       generate.  On  the other hand it  is  not too  efficient  on */
/*       vector  computers   for the  iterative  solution  of  linear */
/*       systems.  Hence, SLAP  changes this input  data structure to */
/*       the SLAP   Column  format for the  iteration (but   does not */
/*       change it back). */

/*       Here is an example of the  SLAP Triad   storage format for a */
/*       5x5 Matrix.  Recall that the entries may appear in any order. */

/*           5x5 Matrix      SLAP Triad format for 5x5 matrix on left. */
/*                              1  2  3  4  5  6  7  8  9 10 11 */
/*       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21 */
/*       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2 */
/*       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       =================== S L A P Column format ================== */

/*       In  this format   the non-zeros are    stored counting  down */
/*       columns (except  for the diagonal  entry, which must  appear */
/*       first in each "column") and are  stored in the real array A. */
/*       In other words,  for  each column    in the matrix   put the */
/*       diagonal  entry  in A.   Then   put  in the  other  non-zero */
/*       elements going   down the  column (except  the  diagonal) in */
/*       order.  The IA array holds the row index  for each non-zero. */
/*       The JA array holds the offsets into the IA, A arrays for the */
/*       beginning   of   each  column.      That is,   IA(JA(ICOL)), */
/*       A(JA(ICOL)) points to the beginning of the ICOL-th column in */
/*       IA and  A.  IA(JA(ICOL+1)-1), A(JA(ICOL+1)-1)  points to the */
/*       end of the ICOL-th column.  Note that we always have JA(N+1) */
/*       = NELT+1, where N is the number of columns in the matrix and */
/*       NELT is the number of non-zeros in the matrix. */

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
/* *See Also: */
/*       SSDCG, SSICCG */
/* ***REFERENCES  1. Louis Hageman \& David Young, ``Applied Iterative */
/*                 Methods'', Academic Press, New York (1981) ISBN */
/*                 0-12-313340-8. */

/*               2. Concus, Golub \& O'Leary, ``A Generalized Conjugate */
/*                 Gradient Method for the Numerical Solution of */
/*                 Elliptic Partial Differential Equations,'' in Sparse */
/*                 Matrix Computations (Bunch \& Rose, Eds.), Academic */
/*                 Press, New York (1979). */
/* ***ROUTINES CALLED  MATVEC, MSOLVE, ISSCG, SCOPY, SDOT, SAXPY, R1MACH */
/* ***END PROLOGUE  SCG */

/*         Check some of the input data. */
/* ***FIRST EXECUTABLE STATEMENT  SCG */
    /* Parameter adjustments */
    --dz;
    --p;
    --z__;
    --r__;
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    *iter = 0;
    *ierr = 0;
    if (*n < 1) {
	*ierr = 3;
	return 0;
    }
    tolmin = r1mach_(&c__3) * 500.f;
    if (*tol < tolmin) {
	*tol = tolmin;
	*ierr = 4;
    }

/*         Calculate initial residual and pseudo-residual, and check */
/*         stopping criterion. */
    (*matvec)(n, &x[1], &r__[1], nelt, &ia[1], &ja[1], &a[1], isym);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__] = b[i__] - r__[i__];
/* L10: */
    }
    (*msolve)(n, &r__[1], &z__[1], nelt, &ia[1], &ja[1], &a[1], isym, &rwork[
	    1], &iwork[1]);

    if (isscg_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)
	    msolve, itol, tol, itmax, iter, err, ierr, iunit, &r__[1], &z__[1]
	    , &p[1], &dz[1], &rwork[1], &iwork[1], &ak, &bk, &bnrm, &solnrm) 
	    != 0) {
	goto L200;
    }
    if (*ierr != 0) {
	return 0;
    }

/*         ***** Iteration loop ***** */

    i__1 = *itmax;
    for (k = 1; k <= i__1; ++k) {
	*iter = k;

/*         Calculate coefficient bk and direction vector p. */
	bknum = sdot_(n, &z__[1], &c__1, &r__[1], &c__1);
	if (bknum <= 0.f) {
	    *ierr = 5;
	    return 0;
	}
	if (*iter == 1) {
	    scopy_(n, &z__[1], &c__1, &p[1], &c__1);
	} else {
	    bk = bknum / bkden;
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		p[i__] = z__[i__] + bk * p[i__];
/* L20: */
	    }
	}
	bkden = bknum;

/*         Calculate coefficient ak, new iterate x, new residual r, */
/*         and new pseudo-residual z. */
	(*matvec)(n, &p[1], &z__[1], nelt, &ia[1], &ja[1], &a[1], isym);
	akden = sdot_(n, &p[1], &c__1, &z__[1], &c__1);
	if (akden <= 0.f) {
	    *ierr = 6;
	    return 0;
	}
	ak = bknum / akden;
	saxpy_(n, &ak, &p[1], &c__1, &x[1], &c__1);
	r__1 = -ak;
	saxpy_(n, &r__1, &z__[1], &c__1, &r__[1], &c__1);
	(*msolve)(n, &r__[1], &z__[1], nelt, &ia[1], &ja[1], &a[1], isym, &
		rwork[1], &iwork[1]);

/*         check stopping criterion. */
	if (isscg_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)
		msolve, itol, tol, itmax, iter, err, ierr, iunit, &r__[1], &
		z__[1], &p[1], &dz[1], &rwork[1], &iwork[1], &ak, &bk, &bnrm, 
		&solnrm) != 0) {
	    goto L200;
	}

/* L100: */
    }

/*         *****   end of loop  ***** */

/*         stopping criterion not satisfied. */
    *iter = *itmax + 1;
    *ierr = 2;

L200:
    return 0;
} /* scg_ */

/* DECK SSDCG */
/* Subroutine */ int ssdcg_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, integer *itol, real 
	*tol, integer *itmax, integer *iter, real *err, integer *ierr, 
	integer *iunit, real *rwork, integer *lenw, integer *iwork, integer *
	leniw)
{
    extern /* Subroutine */ int scg_(integer *, real *, real *, integer *, 
	    integer *, integer *, real *, integer *, S_fp, S_fp, integer *, 
	    real *, integer *, integer *, real *, integer *, integer *, real *
	    , real *, real *, real *, real *, integer *), ss2y_(integer *, 
	    integer *, integer *, integer *, real *, integer *);
    static integer locd, locp, locr;
    extern /* Subroutine */ int ssdi_();
    static integer locw, locz;
    extern /* Subroutine */ int ssds_(integer *, integer *, integer *, 
	    integer *, real *, integer *, real *), ssmv_();
    static integer locdz, lociw;
    extern /* Subroutine */ int schkw_(char *, integer *, integer *, integer *
	    , integer *, integer *, integer *, real *, ftnlen);

/* ***BEGIN PROLOGUE  SSDCG */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSDCG-S), */
/*             Symmetric Linear system, Sparse, Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Diagonally Scaled Conjugate Gradient Sparse Ax=b Solver. */
/*            Routine to solve a  symmetric positive definite  linear */
/*            system Ax  =  b  using   the  Preconditioned  Conjugate */
/*            Gradient   method.    The preconditioner  is   diagonal */
/*            scaling. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX */
/*     INTEGER ITER, IERR, IUNIT, LENW, IWORK(10), LENIW */
/*     REAL B(N), X(N), A(NELT), TOL, ERR, RWORK(5*N) */

/*     CALL SSDCG(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL, */
/*    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW ) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right-hand side vector. */
/* X      :INOUT    Real X(N). */
/*         On input X is your initial guess for solution vector. */
/*         On output X is the final approximate solution. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :INOUT    Integer IA(NELT). */
/* JA     :INOUT    Integer JA(NELT). */
/* A      :INOUT    Integer A(NELT). */
/*         These arrays should hold the matrix A in either the SLAP */
/*         Triad format or the SLAP Column format.  See ``Description'', */
/*         below.  If the SLAP Triad format is chosen it is changed */
/*         internally to the SLAP Column format. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* ITOL   :IN       Integer. */
/*         Flag to indicate type of convergence criterion. */
/*         If ITOL=1, iteration stops when the 2-norm of the residual */
/*         divided by the 2-norm of the right-hand side is less than TOL. */
/*         If ITOL=2, iteration stops when the 2-norm of M-inv times the */
/*         residual divided by the 2-norm of M-inv times the right hand */
/*         side is less than TOL, where M-inv is the inverse of the */
/*         diagonal of A. */
/*         ITOL=11 is often useful for checking and comparing different */
/*         routines.  For this case, the user must supply the "exact" */
/*         solution or a very accurate approximation (one with an error */
/*         much less than TOL) through a common block, */
/*                     COMMON /SOLBLK/ SOLN(1) */
/*         if ITOL=11, iteration stops when the 2-norm of the difference */
/*         between the iterative approximation and the user-supplied */
/*         solution divided by the 2-norm of the user-supplied solution */
/*         is less than TOL.  Note that this requires the user to set up */
/*         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine. */
/*         The routine with this declaration should be loaded before the */
/*         stop test so that the correct length is used by the loader. */
/*         This procedure is not standard Fortran and may not work */
/*         correctly on your system (although it has worked on every */
/*         system the authors have tried).  If ITOL is not 11 then this */
/*         common block is indeed standard Fortran. */
/* TOL    :IN       Real. */
/*         Convergence criterion, as described above. */
/* ITMAX  :IN       Integer. */
/*         Maximum number of iterations. */
/* ITER   :OUT      Integer. */
/*         Number of iterations required to reach convergence, or */
/*         ITMAX+1 if convergence criterion could not be achieved in */
/*         ITMAX iterations. */
/* ERR    :OUT      Real. */
/*         Error estimate of error in final approximate solution, as */
/*         defined by ITOL. */
/* IERR   :OUT      Integer. */
/*         Return error flag. */
/*           IERR = 0 => All went well. */
/*           IERR = 1 => Insufficient storage allocated */
/*                       for WORK or IWORK. */
/*           IERR = 2 => Method failed to converge in */
/*                       ITMAX steps. */
/*           IERR = 3 => Error in user input.  Check input */
/*                       value of N, ITOL. */
/*           IERR = 4 => User error tolerance set too tight. */
/*                       Reset to 500.0*R1MACH(3).  Iteration proceeded. */
/*           IERR = 5 => Preconditioning matrix, M,  is not */
/*                       Positive Definite.  $(r,z) < 0.0$. */
/*           IERR = 6 => Matrix A is not Positive Definite. */
/*                       $(p,Ap) < 0.0$. */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* RWORK  :WORK     Real RWORK(LENW). */
/*         Real array used for workspace. */
/* LENW   :IN       Integer. */
/*         Length of the real workspace, RWORK.  LENW >= 5*N. */
/* IWORK  :WORK     Integer IWORK(LENIW). */
/*         Used to hold pointers into the real workspace, RWORK. */
/*         Upon return the following locations of IWORK hold information */
/*         which may be of use to the user: */
/*         IWORK(9)  Amount of Integer workspace actually used. */
/*         IWORK(10) Amount of Real    workspace actually used. */
/* LENIW  :IN       Integer. */
/*         Length of the integer workspace, IWORK.  LENIW >= 10. */

/* *Description: */
/*       This  routine   performs preconditioned conjugate   gradient */
/*       method on  the  symmetric positive definite   linear  system */
/*       Ax=b.   The preconditioner is  M = DIAG(A), the  diagonal of */
/*       the matrix A.  This is the  simplest of preconditioners  and */
/*       vectorizes very well.   This routine is  simply a driver for */
/*       the SCG routine.  It  calls the SSDS  routine to  set up the */
/*       preconditioning  and  then  calls  SCG  with the appropriate */
/*       MATVEC and MSOLVE routines. */

/*       The Sparse Linear Algebra Package (SLAP) utilizes two matrix */
/*       data structures: 1) the  SLAP Triad  format or  2)  the SLAP */
/*       Column format.  The user can hand this routine either of the */
/*       of these data structures and SLAP  will figure out  which on */
/*       is being used and act accordingly. */

/*       =================== S L A P Triad format =================== */

/*       This routine requires that the  matrix A be   stored in  the */
/*       SLAP  Triad format.  In  this format only the non-zeros  are */
/*       stored.  They may appear in  *ANY* order.  The user supplies */
/*       three arrays of  length NELT, where  NELT is  the number  of */
/*       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For */
/*       each non-zero the user puts the row and column index of that */
/*       matrix element  in the IA and  JA arrays.  The  value of the */
/*       non-zero   matrix  element is  placed  in  the corresponding */
/*       location of the A array.   This is  an  extremely  easy data */
/*       structure to generate.  On  the  other hand it   is  not too */
/*       efficient on vector computers for  the iterative solution of */
/*       linear systems.  Hence,   SLAP changes   this  input    data */
/*       structure to the SLAP Column format  for  the iteration (but */
/*       does not change it back). */

/*       Here is an example of the  SLAP Triad   storage format for a */
/*       5x5 Matrix.  Recall that the entries may appear in any order. */

/*           5x5 Matrix       SLAP Triad format for 5x5 matrix on left. */
/*                              1  2  3  4  5  6  7  8  9 10 11 */
/*       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21 */
/*       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2 */
/*       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       =================== S L A P Column format ================== */

/*       This routine  requires that  the matrix A  be stored in  the */
/*       SLAP Column format.  In this format the non-zeros are stored */
/*       counting down columns (except  for the diagonal entry, which */
/*       must appear first in each  ``column'') and  are stored  in the */
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
/* *Side Effects: */
/*       The SLAP Triad format (IA, JA, A) is modified internally  to */
/*       be the   SLAP   Column format.  See above. */

/* *See Also: */
/*       SCG, SSICCG */
/* ***REFERENCES  1. Louis Hageman \& David Young, ``Applied Iterative */
/*                 Methods'', Academic Press, New York (1981) ISBN */
/*                 0-12-313340-8. */
/*               2. Concus, Golub \& O'Leary, ``A Generalized Conjugate */
/*                 Gradient Method for the Numerical Solution of */
/*                 Elliptic Partial Differential Equations,'' in Sparse */
/*                 Matrix Computations (Bunch \& Rose, Eds.), Academic */
/*                 Press, New York (1979). */
/* ***ROUTINES CALLED  SS2Y, SCHKW, SSDS, SCG */
/* ***END PROLOGUE  SSDCG */

/*         Modify the SLAP matrix data structure to YSMP-Column. */
/* ***FIRST EXECUTABLE STATEMENT  SSDCG */
    /* Parameter adjustments */
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    *ierr = 0;
    if (*n < 1 || *nelt < 1) {
	*ierr = 3;
	return 0;
    }
    ss2y_(n, nelt, &ia[1], &ja[1], &a[1], isym);

/*         Set up the work arrays. */
/*         Compute the inverse of the diagonal of the matrix.  This */
/*         will be used as the preconditioner. */
    lociw = 11;

    locd = 1;
    locr = locd + *n;
    locz = locr + *n;
    locp = locz + *n;
    locdz = locp + *n;
    locw = locdz + *n;

/*         Check the workspace allocations. */
    schkw_("SSDCG", &lociw, leniw, &locw, lenw, ierr, iter, err, (ftnlen)5);
    if (*ierr != 0) {
	return 0;
    }

    iwork[4] = locd;
    iwork[9] = lociw;
    iwork[10] = locw;

    ssds_(n, nelt, &ia[1], &ja[1], &a[1], isym, &rwork[locd]);

/*         Do the Preconditioned Conjugate Gradient. */
    scg_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)ssmv_, (
	    S_fp)ssdi_, itol, tol, itmax, iter, err, ierr, iunit, &rwork[locr]
	    , &rwork[locz], &rwork[locp], &rwork[locdz], &rwork[1], &iwork[1])
	    ;
    return 0;
/* ------------- LAST LINE OF SSDCG FOLLOWS ----------------------------- */
} /* ssdcg_ */

/* DECK SSICCG */
/* Subroutine */ int ssiccg_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, integer *itol, real 
	*tol, integer *itmax, integer *iter, real *err, integer *ierr, 
	integer *iunit, real *rwork, integer *lenw, integer *iwork, integer *
	leniw)
{
    extern /* Subroutine */ int scg_(integer *, real *, real *, integer *, 
	    integer *, integer *, real *, integer *, S_fp, S_fp, integer *, 
	    real *, integer *, integer *, real *, integer *, integer *, real *
	    , real *, real *, real *, real *, integer *);
    static integer nel;
    extern /* Subroutine */ int ss2y_(integer *, integer *, integer *, 
	    integer *, real *, integer *);
    static integer locp, locr, locw, locz;
    extern /* Subroutine */ int ssmv_();
    static integer locel, locdz, lociw;
    extern /* Subroutine */ int schkw_(char *, integer *, integer *, integer *
	    , integer *, integer *, integer *, real *, ftnlen), ssics_(
	    integer *, integer *, integer *, integer *, real *, integer *, 
	    integer *, integer *, integer *, real *, real *, real *, integer *
	    );
    static integer lociel, locjel, locdin;
    extern /* Subroutine */ int ssllti_(), xerrwv_(char *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, integer *, real *, 
	    real *, ftnlen);

/* ***BEGIN PROLOGUE  SSICCG */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSICCG-S), */
/*             Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Incomplete Cholesky */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Incomplete Cholesky Conjugate Gradient Sparse Ax=b Solver. */
/*            Routine to  solve a symmetric  positive definite linear */
/*            system   Ax    =  b  using  the    incomplete  Cholesky */
/*            Preconditioned Conjugate Gradient method. */

/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX */
/*     INTEGER ITER, IERR, IUNIT, LENW, IWORK(NEL+2*n+1), LENIW */
/*     REAL B(N), X(N), A(NELT), TOL, ERR, RWORK(NEL+5*N) */

/*     CALL SSICCG(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL, */
/*    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW ) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right-hand side vector. */
/* X      :INOUT    Real X(N). */
/*         On input X is your initial guess for solution vector. */
/*         On output X is the final approximate solution. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :INOUT    Integer IA(NELT). */
/* JA     :INOUT    Integer JA(NELT). */
/* A      :INOUT    Integer A(NELT). */
/*         These arrays should hold the matrix A in either the SLAP */
/*         Triad format or the SLAP Column format.  See ``Description'', */
/*         below.  If the SLAP Triad format is chosen it is changed */
/*         internally to the SLAP Column format. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* ITOL   :IN       Integer. */
/*         Flag to indicate type of convergence criterion. */
/*         If ITOL=1, iteration stops when the 2-norm of the residual */
/*         divided by the 2-norm of the right-hand side is less than TOL. */
/*         If ITOL=2, iteration stops when the 2-norm of M-inv times the */
/*         residual divided by the 2-norm of M-inv times the right hand */
/*         side is less than TOL, where M-inv is the inverse of the */
/*         diagonal of A. */
/*         ITOL=11 is often useful for checking and comparing different */
/*         routines.  For this case, the user must supply the "exact" */
/*         solution or a very accurate approximation (one with an error */
/*         much less than TOL) through a common block, */
/*                     COMMON /SOLBLK/ SOLN(1) */
/*         if ITOL=11, iteration stops when the 2-norm of the difference */
/*         between the iterative approximation and the user-supplied */
/*         solution divided by the 2-norm of the user-supplied solution */
/*         is less than TOL.  Note that this requires the user to set up */
/*         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine. */
/*         The routine with this declaration should be loaded before the */
/*         stop test so that the correct length is used by the loader. */
/*         This procedure is not standard Fortran and may not work */
/*         correctly on your system (although it has worked on every */
/*         system the authors have tried).  If ITOL is not 11 then this */
/*         common block is indeed standard Fortran. */
/* TOL    :IN       Real. */
/*         Convergence criterion, as described above. */
/* ITMAX  :IN       Integer. */
/*         Maximum number of iterations. */
/* ITER   :OUT      Integer. */
/*         Number of iterations required to reach convergence, or */
/*         ITMAX+1 if convergence criterion could not be achieved in */
/*         ITMAX iterations. */
/* ERR    :OUT      Real. */
/*         Error estimate of error in final approximate solution, as */
/*         defined by ITOL. */
/* IERR   :OUT      Integer. */
/*         Return error flag. */
/*           IERR = 0 => All went well. */
/*           IERR = 1 => Insufficient storage allocated */
/*                       for WORK or IWORK. */
/*           IERR = 2 => Method failed to converge in */
/*                       ITMAX steps. */
/*           IERR = 3 => Error in user input.  Check input */
/*                       value of N, ITOL. */
/*           IERR = 4 => User error tolerance set too tight. */
/*                       Reset to 500.0*R1MACH(3).  Iteration proceeded. */
/*           IERR = 5 => Preconditioning matrix, M,  is not */
/*                       Positive Definite.  $(r,z) < 0.0$. */
/*           IERR = 6 => Matrix A is not Positive Definite. */
/*                       $(p,Ap) < 0.0$. */
/*           IERR = 7 => Incomplete factorization broke down */
/*                       and was fudged.  Resulting preconditioning may */
/*                       be less than the best. */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* RWORK  :WORK     Real RWORK(LENW). */
/*         Real array used for workspace.  NEL is the number of non- */
/*         zeros in the lower triangle of the matrix (including the */
/*         diagonal) */
/* LENW   :IN       Integer. */
/*         Length of the real workspace, RWORK.  LENW >= NEL+5*N. */
/* IWORK  :WORK     Integer IWORK(LENIW). */
/*         Integer array used for workspace.  NEL is the number of non- */
/*         zeros in the lower triangle of the matrix (including the */
/*         diagonal). */
/*         Upon return the following locations of IWORK hold information */
/*         which may be of use to the user: */
/*         IWORK(9)  Amount of Integer workspace actually used. */
/*         IWORK(10) Amount of Real    workspace actually used. */
/* LENIW  :IN       Integer. */
/*         Length of the integer workspace, IWORK.  LENIW >= NEL+N+11. */

/* *Description: */
/*       This routine  performs  preconditioned  conjugate   gradient */
/*       method on the   symmetric positive  definite  linear  system */
/*       Ax=b.   The preconditioner  is  the incomplete Cholesky (IC) */
/*       factorization of the matrix A.  See  SSICS for details about */
/*       the incomplete   factorization algorithm.  One   should note */
/*       here however, that the  IC factorization is a  slow  process */
/*       and  that  one should   save  factorizations  for  reuse, if */
/*       possible.  The   MSOLVE operation (handled  in  SSLLTI) does */
/*       vectorize on machines  with  hardware  gather/scatter and is */
/*       quite fast. */

/*       The Sparse Linear Algebra Package (SLAP) utilizes two matrix */
/*       data structures: 1) the  SLAP Triad  format or  2)  the SLAP */
/*       Column format.  The user can hand this routine either of the */
/*       of these data structures and SLAP  will figure out  which on */
/*       is being used and act accordingly. */

/*       =================== S L A P Triad format =================== */

/*       This routine requires that the  matrix A be   stored in  the */
/*       SLAP  Triad format.  In  this format only the non-zeros  are */
/*       stored.  They may appear in  *ANY* order.  The user supplies */
/*       three arrays of  length NELT, where  NELT is  the number  of */
/*       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For */
/*       each non-zero the user puts the row and column index of that */
/*       matrix element  in the IA and  JA arrays.  The  value of the */
/*       non-zero   matrix  element is  placed  in  the corresponding */
/*       location of the A array.   This is  an  extremely  easy data */
/*       structure to generate.  On  the  other hand it   is  not too */
/*       efficient on vector computers for  the iterative solution of */
/*       linear systems.  Hence,   SLAP changes   this  input    data */
/*       structure to the SLAP Column format  for  the iteration (but */
/*       does not change it back). */

/*       Here is an example of the  SLAP Triad   storage format for a */
/*       5x5 Matrix.  Recall that the entries may appear in any order. */

/*           5x5 Matrix       SLAP Triad format for 5x5 matrix on left. */
/*                              1  2  3  4  5  6  7  8  9 10 11 */
/*       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21 */
/*       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2 */
/*       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1 */
/*       | 0  0  0 44  0| */
/*       |51  0 53  0 55| */

/*       =================== S L A P Column format ================== */

/*       This routine  requires that  the matrix A  be stored in  the */
/*       SLAP Column format.  In this format the non-zeros are stored */
/*       counting down columns (except  for the diagonal entry, which */
/*       must appear first in each ``column'') and  are stored  in the */
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
/* *Side Effects: */
/*       The SLAP Triad format (IA, JA, A) is modified internally to be */
/*       the SLAP Column format.  See above. */

/* *See Also: */
/*       SCG, SSLLTI */
/* ***REFERENCES  1. Louis Hageman \& David Young, ``Applied Iterative */
/*                 Methods'', Academic Press, New York (1981) ISBN */
/*                 0-12-313340-8. */
/*               2. Concus, Golub \& O'Leary, ``A Generalized Conjugate */
/*                 Gradient Method for the Numerical Solution of */
/*                 Elliptic Partial Differential Equations,'' in Sparse */
/*                 Matrix Computations (Bunch \& Rose, Eds.), Academic */
/*                 Press, New York (1979). */
/* ***ROUTINES CALLED  SS2Y, SCHKW, SSICS, XERRWV, SCG */
/* ***END PROLOGUE  SSICCG */

/*         Change the SLAP input matrix IA, JA, A to SLAP-Column format. */
/* ***FIRST EXECUTABLE STATEMENT  SSICCG */
    /* Parameter adjustments */
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    *ierr = 0;
    if (*n < 1 || *nelt < 1) {
	*ierr = 3;
	return 0;
    }
    ss2y_(n, nelt, &ia[1], &ja[1], &a[1], isym);

/*         Count number of elements in lower triangle of the matrix. */
/*         Then set up the work arrays. */
    if (*isym == 0) {
	nel = (*nelt + *n) / 2;
    } else {
	nel = *nelt;
    }

    locjel = 11;
    lociel = locjel + nel;
    lociw = lociel + *n + 1;

    locel = 1;
    locdin = locel + nel;
    locr = locdin + *n;
    locz = locr + *n;
    locp = locz + *n;
    locdz = locp + *n;
    locw = locdz + *n;

/*         Check the workspace allocations. */
    schkw_("SSICCG", &lociw, leniw, &locw, lenw, ierr, iter, err, (ftnlen)6);
    if (*ierr != 0) {
	return 0;
    }

    iwork[1] = nel;
    iwork[2] = locjel;
    iwork[3] = lociel;
    iwork[4] = locel;
    iwork[5] = locdin;
    iwork[9] = lociw;
    iwork[10] = locw;

/*         Compute the Incomplete Cholesky decomposition. */

    ssics_(n, nelt, &ia[1], &ja[1], &a[1], isym, &nel, &iwork[lociel], &iwork[
	    locjel], &rwork[locel], &rwork[locdin], &rwork[locr], ierr);
    if (*ierr != 0) {
	xerrwv_("SSICCG: Warning...IC factorization broke down on step i1.  "
		"Diagonal was set to unity and factorization proceeded.", &
		c__113, &c__1, &c__1, &c__1, ierr, &c__0, &c__0, &c_b28, &
		c_b28, (ftnlen)113);
	*ierr = 7;
    }

/*         Do the Preconditioned Conjugate Gradient. */
    scg_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)ssmv_, (
	    S_fp)ssllti_, itol, tol, itmax, iter, err, ierr, iunit, &rwork[
	    locr], &rwork[locz], &rwork[locp], &rwork[locdz], &rwork[1], &
	    iwork[1]);
    return 0;
/* ------------- LAST LINE OF SSICCG FOLLOWS ---------------------------- */
} /* ssiccg_ */

/* DECK ISSCG */
integer isscg_(integer *n, real *b, real *x, integer *nelt, integer *ia, 
	integer *ja, real *a, integer *isym, S_fp msolve, integer *itol, real 
	*tol, integer *itmax, integer *iter, real *err, integer *ierr, 
	integer *iunit, real *r__, real *z__, real *p, real *dz, real *rwork, 
	integer *iwork, real *ak, real *bk, real *bnrm, real *solnrm)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 Preconditioned Conjugate Gradient for"
	    " \002,\002N, ITOL = \002,i5,i5,/\002 ITER\002,\002   Error Estim"
	    "ate\002,\002            Alpha\002,\002             Beta\002)";
    static char fmt_1010[] = "(1x,i4,1x,e16.7,1x,e16.7,1x,e16.7)";

    /* System generated locals */
    integer ret_val, i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__;
    extern doublereal snrm2_(integer *, real *, integer *);

    /* Fortran I/O blocks */
    static cilist io___30 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_1010, 0 };


/* ***BEGIN PROLOGUE  ISSCG */
/* ***REFER TO  SCG, SSDCG, SSICCG */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(ISSCG-S), */
/*             Linear system, Sparse, Stop Test */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Preconditioned Conjugate Gradient Stop Test. */
/*            This routine calculates the stop test for the Conjugate */
/*            Gradient iteration scheme.  It returns a nonzero if the */
/*            error estimate (the type of which is determined by ITOL) */
/*            is less than the user specified tolerance TOL. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER */
/*     INTEGER IERR, IUNIT, IWORK(USER DEFINED) */
/*     REAL    B(N), X(N), A(N), TOL, ERR, R(N), Z(N), P(N), DZ(N) */
/*     REAL    RWORK(USER DEFINED), AK, BK, BNRM, SOLNRM */
/*     EXTERNAL MSOLVE */

/*     IF( ISSCG(N, B, X, NELT, IA, JA, A, ISYM, MSOLVE, ITOL, TOL, */
/*    $     ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, DZ, RWORK, IWORK, */
/*    $     AK, BK, BNRM, SOLNRM) .NE. 0 ) THEN ITERATION DONE */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right-hand side vector. */
/* X      :IN       Real X(N). */
/*         The current approximate solution vector. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Real A(NELT). */
/*         These arrays should hold the matrix A in either the SLAP */
/*         Triad format or the SLAP Column format.  See ``Description'' */
/*         in the SCG, SSDCG or SSICCG routines. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* MSOLVE :EXT      External. */
/*         Name of a routine which solves a linear system MZ = R for */
/*         Z given R with the preconditioning matrix M (M is supplied via */
/*         RWORK and IWORK arrays).  The name of the MSOLVE routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MSLOVE is: */
/*             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK) */
/*         Where N is the number of unknowns, R is the right-hand side */
/*         vector, and Z is the solution upon return.  RWORK is a real */
/*         array that can be used to pass necessary preconditioning */
/*         information and/or workspace to MSOLVE.  IWORK is an integer */
/*         work array for the same purpose as RWORK. */
/* ITOL   :IN       Integer. */
/*         Flag to indicate type of convergence criterion. */
/*         If ITOL=1, iteration stops when the 2-norm of the residual */
/*         divided by the 2-norm of the right-hand side is less than TOL. */
/*         If ITOL=2, iteration stops when the 2-norm of M-inv times the */
/*         residual divided by the 2-norm of M-inv times the right hand */
/*         side is less than tol, where M-inv is the inverse of the */
/*         diagonal of A. */
/*         ITOL=11 is often useful for checking and comparing different */
/*         routines.  For this case, the user must supply the ``exact'' */
/*         solution or a very accurate approximation (one with an error */
/*         much less than tol) through a common block, */
/*         COMMON /SOLBLK/ SOLN( ) */
/*         if ITOL=11, iteration stops when the 2-norm of the difference */
/*         between the iterative approximation and the user-supplied */
/*         solution divided by the 2-norm of the user-supplied solution */
/*         is less than tol. */
/* TOL    :IN       Real. */
/*         Convergence criterion, as described above. */
/* ITMAX  :IN       Integer. */
/*         Maximum number of iterations. */
/* ITER   :IN       Integer. */
/*         The iteration for which to check for convergence. */
/* ERR    :OUT      Real. */
/*         Error estimate of error in the X(N) approximate solution, as */
/*         defined by ITOL. */
/* IERR   :OUT      Integer. */
/*         Error flag.  IERR is set to 3 if ITOL is not on of the */
/*         acceptable values, see above. */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* R      :IN       Real R(N). */
/*         The residual R = B-AX. */
/* Z      :WORK     Real Z(N). */
/*         Workspace used to hold the pseudo-residual M Z = R. */
/* P      :IN       Real P(N). */
/*         The conjugate direction vector. */
/* DZ     :WORK     Real DZ(N). */
/*         Workspace used to hold temporary vector(s). */
/* RWORK  :WORK     Real RWORK(USER DEFINABLE). */
/*         Real array that can be used by MSOLVE. */
/* IWORK  :WORK     Integer IWORK(USER DEFINABLE). */
/*         Integer array that can be used by MSOLVE. */
/* BNRM   :INOUT    Real. */
/*         Norm of the right hand side.  Type of norm depends on ITOL. */
/*         Calculated only on the first call. */
/* SOLNRM :INOUT    Real. */
/*         2-Norm of the true solution, SOLN.  Only computed and used */
/*         if ITOL = 11. */

/* *Function Return Values: */
/*       0 : Error estimate (determined by ITOL) is *NOT* less than the */
/*           specified tolerance, TOL.  The iteration must continue. */
/*       1 : Error estimate (determined by ITOL) is less than the */
/*           specified tolerance, TOL.  The iteration can be considered */
/*           complete. */

/* *Precision:           Single Precision */
/* *See Also: */
/*       SCG, SSDCG, SSICCG */

/* *Cautions: */
/*     This routine will attempt to write to the fortran logical output */
/*     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that */
/*     this  logical  unit  must  be  attached  to  a  file or terminal */
/*     before calling this routine with a non-zero  value  for   IUNIT. */
/*     This routine does not check for the validity of a non-zero IUNIT */
/*     unit number. */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  MSOLVE, SNRM2 */
/* ***COMMON BLOCKS    SOLBLK */
/* ***END PROLOGUE  ISSCG */

/* ***FIRST EXECUTABLE STATEMENT  ISSCG */
    /* Parameter adjustments */
    --dz;
    --p;
    --z__;
    --r__;
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rwork;
    --iwork;

    /* Function Body */
    ret_val = 0;

    if (*itol == 1) {
/*         err = ||Residual||/||RightHandSide|| (2-Norms). */
	if (*iter == 0) {
	    *bnrm = snrm2_(n, &b[1], &c__1);
	}
	*err = snrm2_(n, &r__[1], &c__1) / *bnrm;
    } else if (*itol == 2) {
/*                  -1              -1 */
/*         err = ||M  Residual||/||M  RightHandSide|| (2-Norms). */
	if (*iter == 0) {
	    (*msolve)(n, &r__[1], &z__[1], nelt, &ia[1], &ja[1], &a[1], isym, 
		    &rwork[1], &iwork[1]);
	    *bnrm = snrm2_(n, &dz[1], &c__1);
	}
	*err = snrm2_(n, &z__[1], &c__1) / *bnrm;
    } else if (*itol == 11) {
/*         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms). */
	if (*iter == 0) {
	    *solnrm = snrm2_(n, solblk_1.soln, &c__1);
	}
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dz[i__] = x[i__] - solblk_1.soln[i__ - 1];
/* L10: */
	}
	*err = snrm2_(n, &dz[1], &c__1) / *solnrm;
    } else {

/*         If we get here ITOL is not one of the acceptable values. */
	*err = 1e10f;
	*ierr = 3;
    }

    if (*iunit != 0) {
	if (*iter == 0) {
	    io___30.ciunit = *iunit;
	    s_wsfe(&io___30);
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*itol), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	io___31.ciunit = *iunit;
	s_wsfe(&io___31);
	do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*err), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*ak), (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*bk), (ftnlen)sizeof(real));
	e_wsfe();
    }
    if (*err <= *tol) {
	ret_val = 1;
    }
    return ret_val;
/* ------------- LAST LINE OF ISSCG FOLLOWS ------------------------------ */
} /* isscg_ */

