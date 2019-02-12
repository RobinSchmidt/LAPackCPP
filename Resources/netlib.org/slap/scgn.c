/* scgn.f -- translated by f2c (version 20100827).
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

/* DECK SCGN */
/* Subroutine */ int scgn_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, S_fp matvec, S_fp 
	mttvec, S_fp msolve, integer *itol, real *tol, integer *itmax, 
	integer *iter, real *err, integer *ierr, integer *iunit, real *r__, 
	real *z__, real *p, real *atp, real *atz, real *dz, real *atdz, real *
	rwork, integer *iwork)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Local variables */
    static integer i__, k;
    static real ak, bk, bnrm;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real akden, bkden, bknum;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), saxpy_(integer *, real *, real *, integer *, real *, 
	    integer *);
    extern doublereal r1mach_(integer *);
    extern integer isscgn_(integer *, real *, real *, integer *, integer *, 
	    integer *, real *, integer *, S_fp, S_fp, S_fp, integer *, real *,
	     integer *, integer *, real *, integer *, integer *, real *, real 
	    *, real *, real *, real *, real *, real *, real *, integer *, 
	    real *, real *, real *, real *);
    static real tolmin, solnrm;

/* ***BEGIN PROLOGUE  SCGN */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SCGN-S), */
/*             Non-Symmetric Linear system solve, Sparse, */
/*             Iterative Precondition, Normal Equations. */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Preconditioned CG Sparse Ax=b Solver for Normal Equations. */
/*            Routine  to solve a general linear system Ax = b using the */
/*            Preconditioned Conjugate Gradient method  applied to   the */
/*            normal equations AA'y = b, x=A'y. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX */
/*     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINABLE) */
/*     REAL     B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N), P(N), ATP(N) */
/*     REAL     ATZ(N), DZ(N), ATDZ(N) */
/*     REAL     RWORK(USER DEFINABLE) */
/*     EXTERNAL MATVEC, MTTVEC, MSOLVE */

/*     CALL SCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, */
/*    $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, */
/*    $     Z, P, ATP, ATZ, DZ, ATDZ, RWORK, IWORK) */

/* *Arguments: */
/* N      :IN       Integer */
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
/*         It could take any form.  See "Description", below */
/*         for more late breaking details... */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* MATVEC :EXT      External. */
/*         Name of a routine which performs the matrix vector multiply */
/*         y = A*X given A and X.  The name of the MATVEC routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MATVEC is: */
/*             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM ) */
/*         Where N is the number of unknowns, Y is the product A*X */
/*         upon return X is an input vector, NELT is the number of */
/*         non-zeros in the SLAP-Column IA, JA, A storage for the matrix */
/*         A.  ISYM is a flag which, if non-zero, denotes that A is */
/*         symmetric and only the lower or upper triangle is stored. */
/* MTTVEC :EXT      External. */
/*         Name of a routine which performs the matrix transpose vector */
/*         multiply y = A'*X given A and X (where ' denotes transpose). */
/*         The name of the MTTVEC routine must be declared external in */
/*         the calling program.  The calling sequence to MTTVEC is the */
/*         same as that for MATVEC, viz.: */
/*             CALL MTTVEC( N, X, Y, NELT, IA, JA, A, ISYM ) */
/*         Where N is the number of unknowns, Y is the product A'*X */
/*         upon return X is an input vector, NELT is the number of */
/*         non-zeros in the SLAP-Column IA, JA, A storage for the matrix */
/*         A.  ISYM is a flag which, if non-zero, denotes that A is */
/*         symmetric and only the lower or upper triangle is stored. */
/* MSOLVE :EXT      External. */
/*         Name of a routine which solves a linear system MZ = R for */
/*         Z given R with the preconditioning matrix M (M is supplied via */
/*         RWORK and IWORK arrays).  The name of the MSOLVE routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MSOLVE is: */
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
/* ATP    :WORK     Real ATP(N). */
/* ATZ    :WORK     Real ATZ(N). */
/* DZ     :WORK     Real DZ(N). */
/* ATDZ   :WORK     Real ATDZ(N). */
/* RWORK  :WORK     Real RWORK(USER DEFINABLE). */
/*         Real array that can be used by  MSOLVE. */
/* IWORK  :WORK     Integer IWORK(USER DEFINABLE). */
/*         Integer array that can be used by  MSOLVE. */

/* *Description: */
/*       This  routine applies the  preconditioned conjugate gradient */
/*       (PCG) method to a non-symmetric system of equations Ax=b. To */
/*       do this the normal equations are solved: */
/*               AA' y  = b, where  x  = A'y. */
/*       In PCG method the iteration count is determined by condition */
/*                               -1 */
/*       number of the  matrix (M  A).   In the  situation where  the */
/*       normal equations are  used  to solve a  non-symmetric system */
/*       the condition number depends on  AA' and should therefore be */
/*       much worse than that of A.  This is the conventional wisdom. */
/*       When one has a good preconditioner for AA' this may not hold. */
/*       The latter is the situation when SCGN should be tried. */

/*       If one is trying to solve  a symmetric system, SCG should be */
/*       used instead. */

/*       This routine does  not care  what matrix data   structure is */
/*       used for  A and M.  It simply   calls  the MATVEC and MSOLVE */
/*       routines, with  the arguments as  described above.  The user */
/*       could write any type of structure and the appropriate MATVEC */
/*       and MSOLVE routines.  It is assumed  that A is stored in the */
/*       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is */
/*       stored  in  IWORK  and  RWORK)  in  some fashion.   The SLAP */
/*       routines SSDCGN and SSLUCN are examples of this procedure. */

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

/*           5x5 Matrix       SLAP Triad format for 5x5 matrix on left. */
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
/*         SSDCGN, SSLUCN, ISSCGN */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  MATVEC, MTTVEC, MSOLVE, ISSCGN, */
/*                    SCOPY, SDOT, SAXPY, R1MACH */
/* ***END PROLOGUE  SCGN */

/*         Check user input. */
/* ***FIRST EXECUTABLE STATEMENT  SCGN */
    /* Parameter adjustments */
    --atdz;
    --dz;
    --atz;
    --atp;
    --p;
    --z__;
    --r__;
    --a;
    --x;
    --b;
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
    (*mttvec)(n, &z__[1], &atz[1], nelt, &ia[1], &ja[1], &a[1], isym);

    if (isscgn_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)
	    matvec, (S_fp)mttvec, (S_fp)msolve, itol, tol, itmax, iter, err, 
	    ierr, iunit, &r__[1], &z__[1], &p[1], &atp[1], &atz[1], &dz[1], &
	    atdz[1], &rwork[1], &iwork[1], &ak, &bk, &bnrm, &solnrm) != 0) {
	goto L200;
    }
    if (*ierr != 0) {
	return 0;
    }

/*         ***** iteration loop ***** */

    i__1 = *itmax;
    for (k = 1; k <= i__1; ++k) {
	*iter = k;

/*         Calculate coefficient BK and direction vector P. */
	bknum = sdot_(n, &z__[1], &c__1, &r__[1], &c__1);
	if (bknum <= 0.f) {
	    *ierr = 6;
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

/*         Calculate coefficient AK, new iterate X, new residual R, */
/*         and new pseudo-residual ATZ. */
	if (*iter != 1) {
	    saxpy_(n, &bk, &atp[1], &c__1, &atz[1], &c__1);
	}
	scopy_(n, &atz[1], &c__1, &atp[1], &c__1);
	akden = sdot_(n, &atp[1], &c__1, &atp[1], &c__1);
	if (akden <= 0.f) {
	    *ierr = 6;
	    return 0;
	}
	ak = bknum / akden;
	saxpy_(n, &ak, &atp[1], &c__1, &x[1], &c__1);
	(*matvec)(n, &atp[1], &z__[1], nelt, &ia[1], &ja[1], &a[1], isym);
	r__1 = -ak;
	saxpy_(n, &r__1, &z__[1], &c__1, &r__[1], &c__1);
	(*msolve)(n, &r__[1], &z__[1], nelt, &ia[1], &ja[1], &a[1], isym, &
		rwork[1], &iwork[1]);
	(*mttvec)(n, &z__[1], &atz[1], nelt, &ia[1], &ja[1], &a[1], isym);

/*         check stopping criterion. */
	if (isscgn_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)
		matvec, (S_fp)mttvec, (S_fp)msolve, itol, tol, itmax, iter, 
		err, ierr, iunit, &r__[1], &z__[1], &p[1], &atp[1], &atz[1], &
		dz[1], &atdz[1], &rwork[1], &iwork[1], &ak, &bk, &bnrm, &
		solnrm) != 0) {
	    goto L200;
	}

/* L100: */
    }

/*         *****   end of loop  ***** */

/*         stopping criterion not satisfied. */
    *iter = *itmax + 1;

L200:
    return 0;
/* ------------- LAST LINE OF SCGN FOLLOWS ---------------------------- */
} /* scgn_ */

/* DECK SSDCGN */
/* Subroutine */ int ssdcgn_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, integer *itol, real 
	*tol, integer *itmax, integer *iter, real *err, integer *ierr, 
	integer *iunit, real *rwork, integer *lenw, integer *iwork, integer *
	leniw)
{
    extern /* Subroutine */ int ss2y_(integer *, integer *, integer *, 
	    integer *, real *, integer *);
    static integer locd;
    extern /* Subroutine */ int scgn_(integer *, real *, real *, integer *, 
	    integer *, integer *, real *, integer *, S_fp, S_fp, S_fp, 
	    integer *, real *, integer *, integer *, real *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *, integer *);
    static integer locp, locr;
    extern /* Subroutine */ int ssdi_();
    static integer locw, locz;
    extern /* Subroutine */ int ssmv_(), ssd2s_(integer *, integer *, integer 
	    *, integer *, real *, integer *, real *);
    static integer locdz, lociw;
    extern /* Subroutine */ int schkw_(char *, integer *, integer *, integer *
	    , integer *, integer *, integer *, real *, ftnlen), ssmtv_();
    static integer locatd, locatp, locatz;

/* ***BEGIN PROLOGUE  SSDCGN */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSDCGN-S), */
/*             Non-Symmetric Linear system solve, Sparse, */
/*             Iterative Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Diagonally Scaled CG Sparse Ax=b Solver for Normal Eqn's. */
/*            Routine to solve a general linear system Ax = b using */
/*            diagonal scaling with the Conjugate  Gradient  method */
/*            applied to the the normal equations, viz.,  AA'y = b, */
/*            where x = A'y. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX */
/*     INTEGER ITER, IERR, IUNIT, LENW, IWORK, LENIW */
/*     REAL B(N), X(N), A(NELT), TOL, ERR, RWORK(8*N) */

/*     CALL SSDCGN(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL, */
/*    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW) */

/* *Arguments: */
/* N      :IN       Integer */
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
/*         Triad format or the SLAP Column format.  See "Description", */
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
/*         Length of the real workspace, RWORK.  LENW >= 8*N. */
/* IWORK  :WORK     Integer IWORK(LENIW). */
/*         Used to hold pointers into the RWORK array. */
/*         Upon return the following locations of IWORK hold information */
/*         which may be of use to the user: */
/*         IWORK(9)  Amount of Integer workspace actually used. */
/*         IWORK(10) Amount of Real    workspace actually used. */
/* LENIW  :IN       Integer. */
/*         Length of the integer workspace, IWORK.  LENIW >= 10. */

/* *Description: */
/*       This  routine is simply a driver  for the  SCGN routine.  It */
/*       calls the   SSD2S  routine to set up the preconditioning and */
/*       then calls SCGN with the appropriate   MATVEC  and    MSOLVE */
/*       routines. */

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

/*           5x5 Matrix      SLAP Triad format for 5x5 matrix on left. */
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
/* *Side Effects: */
/*       The SLAP Triad format (IA, JA, A) is modified internally to be */
/*       the SLAP Column format.  See above. */

/* *See Also: */
/*       SCGN, SSD2S, SSMV, SSMTV, SSDI */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  SS2Y, SCHKW, SSD2S, SCGN, SSMV, SSMTV, SSDI */
/* ***END PROLOGUE  SSDCGN */

/*         Modify the SLAP matrix data structure to YSMP-Column. */
/* ***FIRST EXECUTABLE STATEMENT  SSDCGN */
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
/*         Compute the inverse of the diagonal of AA'.  This will be */
/*         used as the preconditioner. */
    lociw = 11;

    locd = 1;
    locr = locd + *n;
    locz = locr + *n;
    locp = locz + *n;
    locatp = locp + *n;
    locatz = locatp + *n;
    locdz = locatz + *n;
    locatd = locdz + *n;
    locw = locatd + *n;

/*         Check the workspace allocations. */
    schkw_("SSDCGN", &lociw, leniw, &locw, lenw, ierr, iter, err, (ftnlen)6);
    if (*ierr != 0) {
	return 0;
    }

    iwork[4] = locd;
    iwork[9] = lociw;
    iwork[10] = locw;

    ssd2s_(n, nelt, &ia[1], &ja[1], &a[1], isym, &rwork[1]);

/*         Perform Conjugate Gradient algorithm on the normal equations. */
    scgn_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)ssmv_, (
	    S_fp)ssmtv_, (S_fp)ssdi_, itol, tol, itmax, iter, err, ierr, 
	    iunit, &rwork[locr], &rwork[locz], &rwork[locp], &rwork[locatp], &
	    rwork[locatz], &rwork[locdz], &rwork[locatd], &rwork[1], &iwork[1]
	    );

    if (*iter > *itmax) {
	*ierr = 2;
    }
    return 0;
/* ------------- LAST LINE OF SSDCGN FOLLOWS ---------------------------- */
} /* ssdcgn_ */

/* DECK SSLUCN */
/* Subroutine */ int sslucn_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, integer *itol, real 
	*tol, integer *itmax, integer *iter, real *err, integer *ierr, 
	integer *iunit, real *rwork, integer *lenw, integer *iwork, integer *
	leniw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, nl, nu;
    extern /* Subroutine */ int ss2y_(integer *, integer *, integer *, 
	    integer *, real *, integer *);
    static integer jbgn, jend, icol, locl;
    extern /* Subroutine */ int scgn_(integer *, real *, real *, integer *, 
	    integer *, integer *, real *, integer *, S_fp, S_fp, S_fp, 
	    integer *, real *, integer *, integer *, real *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *, real *,
	     real *, integer *);
    static integer locp, locr, locu, locw, locz;
    extern /* Subroutine */ int ssmv_();
    static integer locnc, locil, locjl, lociu, locju, locnr, lociw, locdz;
    extern /* Subroutine */ int schkw_(char *, integer *, integer *, integer *
	    , integer *, integer *, integer *, real *, ftnlen), ssmtv_();
    static integer locatd, locdin, locatp, locatz;
    extern /* Subroutine */ int ssmmti_(), ssilus_(integer *, integer *, 
	    integer *, integer *, real *, integer *, integer *, integer *, 
	    integer *, real *, real *, integer *, integer *, integer *, real *
	    , integer *, integer *);

/* ***BEGIN PROLOGUE  SSLUCN */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSLUCN-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Incomplete LU Precondition */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Incomplete LU CG Sparse Ax=b Solver for Normal Equations. */
/*            Routine to solve  a general linear system Ax = b using the */
/*            incomplete  LU decomposition  with  the Conjugate Gradient */
/*            method  applied to the normal equations,  viz., AA'y =  b, */
/*            x=A'y. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX */
/*     INTEGER ITER, IERR, IUNIT, LENW, IWORK(NEL+NU+4*N+2), LENIW */
/*     REAL B(N), X(N), A(NELT), TOL, ERR, RWORK(NEL+NU+8*N) */

/*     CALL SSLUCN(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL, */
/*    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW ) */

/* *Arguments: */
/* N      :IN       Integer */
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
/*         Triad format or the SLAP Column format.  See "Description", */
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
/*         diagonal).  NU is the number of nonzeros in the upper */
/*         triangle of the matrix (including the diagonal). */
/* LENW   :IN       Integer. */
/*         Length of the real workspace, RWORK.  LENW >= NEL+NU+8*N. */
/* IWORK  :WORK     Integer IWORK(LENIW). */
/*         Integer array used for workspace.  NEL is the number of non- */
/*         zeros in the lower triangle of the matrix (including the */
/*         diagonal).  NU is the number of nonzeros in the upper */
/*         triangle of the matrix (including the diagonal). */
/*         Upon return the following locations of IWORK hold information */
/*         which may be of use to the user: */
/*         IWORK(9)  Amount of Integer workspace actually used. */
/*         IWORK(10) Amount of Real    workspace actually used. */
/* LENIW  :IN       Integer. */
/*         Length of the integer workspace, IWORK.  LENIW >= */
/*         NEL+NU+4*N+12. */

/* *Description: */
/*       This  routine is simply a driver  for the  SCGN  routine.    It */
/*       calls the SSILUS routine to set up the preconditioning and then */
/*       calls SCGN with the appropriate  MATVEC  and  MSOLVE  routines. */

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
/* *Side Effects: */
/*       The SLAP Triad format (IA, JA, A) is modified internally to be */
/*       the SLAP Column format.  See above. */

/* *See Also: */
/*       SCGN, SDCGN, SSILUS */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  SS2Y, SSILUS, SCHKW, SSMV, SSMTV, SSMMTI, SCGN */
/* ***END PROLOGUE  SSLUCN */


/*         Change the SLAP input matrix IA, JA, A to SLAP-Column format. */
/* ***FIRST EXECUTABLE STATEMENT  SSLUCN */
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

/*         Count number of Non-Zero elements preconditioner ILU matrix. */
/*         Then set up the work arrays. */
    nl = 0;
    nu = 0;
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
/*         Don't count diagional. */
	jbgn = ja[icol] + 1;
	jend = ja[icol + 1] - 1;
	if (jbgn <= jend) {
/* VD$ NOVECTOR */
	    i__2 = jend;
	    for (j = jbgn; j <= i__2; ++j) {
		if (ia[j] > icol) {
		    ++nl;
		    if (*isym != 0) {
			++nu;
		    }
		} else {
		    ++nu;
		}
/* L10: */
	    }
	}
/* L20: */
    }

    locil = 11;
    locjl = locil + *n + 1;
    lociu = locjl + nl;
    locju = lociu + nu;
    locnr = locju + *n + 1;
    locnc = locnr + *n;
    lociw = locnc + *n;

    locl = 1;
    locdin = locl + nl;
    locu = locdin + *n;
    locr = locu + nu;
    locz = locr + *n;
    locp = locz + *n;
    locatp = locp + *n;
    locatz = locatp + *n;
    locdz = locatz + *n;
    locatd = locdz + *n;
    locw = locatd + *n;

/*         Check the workspace allocations. */
    schkw_("SSLUCN", &lociw, leniw, &locw, lenw, ierr, iter, err, (ftnlen)6);
    if (*ierr != 0) {
	return 0;
    }

    iwork[1] = locil;
    iwork[2] = locjl;
    iwork[3] = lociu;
    iwork[4] = locju;
    iwork[5] = locl;
    iwork[6] = locdin;
    iwork[7] = locu;
    iwork[9] = lociw;
    iwork[10] = locw;

/*         Compute the Incomplete LU decomposition. */
    ssilus_(n, nelt, &ia[1], &ja[1], &a[1], isym, &nl, &iwork[locil], &iwork[
	    locjl], &rwork[locl], &rwork[locdin], &nu, &iwork[lociu], &iwork[
	    locju], &rwork[locu], &iwork[locnr], &iwork[locnc]);

/*         Perform Conjugate Gradient algorithm on the normal equations. */
    scgn_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)ssmv_, (
	    S_fp)ssmtv_, (S_fp)ssmmti_, itol, tol, itmax, iter, err, ierr, 
	    iunit, &rwork[locr], &rwork[locz], &rwork[locp], &rwork[locatp], &
	    rwork[locatz], &rwork[locdz], &rwork[locatd], &rwork[1], &iwork[1]
	    );

    if (*iter > *itmax) {
	*ierr = 2;
    }
    return 0;
/* ------------- LAST LINE OF SSLUCN FOLLOWS ---------------------------- */
} /* sslucn_ */

/* DECK ISSCGN */
integer isscgn_(integer *n, real *b, real *x, integer *nelt, integer *ia, 
	integer *ja, real *a, integer *isym, S_fp matvec, S_fp mttvec, S_fp 
	msolve, integer *itol, real *tol, integer *itmax, integer *iter, real 
	*err, integer *ierr, integer *iunit, real *r__, real *z__, real *p, 
	real *atp, real *atz, real *dz, real *atdz, real *rwork, integer *
	iwork, real *ak, real *bk, real *bnrm, real *solnrm)
{
    /* Format strings */
    static char fmt_1000[] = "(\002 PCG Applied to the Normal Equations for"
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
    static cilist io___46 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___47 = { 0, 0, 0, fmt_1010, 0 };


/* ***BEGIN PROLOGUE  ISSCGN */
/* ***REFER TO  SCGN, SSDCGN, SSLUCN */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(ISSCGN-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Normal Equations */
/* ***AUTHOR  Greenbaum, Anne, Courant Institute */
/*           Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Preconditioned CG on Normal Equations Stop Test. */
/*            This routine calculates the stop test for the Conjugate */
/*            Gradient iteration   scheme applied  to     the  normal */
/*            equations.  It returns a nonzero  if the error estimate */
/*            (the type of which is determined by  ITOL) is less than */
/*            the user specified tolerance TOL. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER */
/*     INTEGER  IERR, IUNIT, IWORK(USER DEFINED) */
/*     REAL     B(N), X(N), A(N), TOL, ERR, R(N), Z(N), P(N), ATP(N) */
/*     REAL     ATZ(N), DZ(N), ATDZ(N) */
/*     REAL     RWORK(USER DEFINED), AK, BK, BNRM, SOLNRM */
/*     EXTERNAL MATVEC, MTTVEC, MSOLVE */

/*     IF( ISTPCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, */
/*    $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, */
/*    $     ATP, ATZ, DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM, SOLNRM) */
/*    $     .NE. 0 ) THEN ITERATION DONE */

/* *Arguments: */
/* N      :IN       Integer */
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
/*         These arrays contain the matrix data structure for A. */
/*         It could take any form.  See "Description" in the */
/*         SSCGN routine. */
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
/*         non-zeros in the SLAP-Column IA, JA, A storage for the matrix */
/*         A.  ISYM is a flag which, if non-zero, denotes that A is */
/*         symmetric and only the lower or upper triangle is stored. */
/* MTTVEC :EXT      External. */
/*         Name of a routine which performs the matrix transpose vector */
/*         multiply y = A'*X given A and X (where ' denotes transpose). */
/*         The name of the MTTVEC routine must be declared external in */
/*         the calling program.  The calling sequence to MTTVEC is the */
/*         same as that for MATVEC, viz.: */
/*             CALL MTTVEC( N, X, Y, NELT, IA, JA, A, ISYM ) */
/*         Where N is the number of unknowns, Y is the product A'*X */
/*         upon return X is an input vector, NELT is the number of */
/*         non-zeros in the SLAP-Column IA, JA, A storage for the matrix */
/*         A.  ISYM is a flag which, if non-zero, denotes that A is */
/*         symmetric and only the lower or upper triangle is stored. */
/* MSOLVE :EXT      External. */
/*         Name of a routine which solves a linear system MZ = R for */
/*         Z given R with the preconditioning matrix M (M is supplied via */
/*         RWORK and IWORK arrays).  The name of the MSOLVE routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MSOLVE is: */
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
/* P      :IN       Real P(N). */
/*         The conjugate direction vector. */
/* ATP    :IN       Real ATP(N). */
/*         A-transpose times the conjugate direction vector. */
/* ATZ    :IN       Real ATZ(N). */
/*         A-transpose times the pseudo-residual. */
/* DZ     :IN       Real DZ(N). */
/*         Workspace used to hold temporary vector(s). */
/* ATDZ   :WORK       Real ATDZ(N). */
/*         Workspace. */
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
/*       SSCGN */

/* *Cautions: */
/*     This routine will attempt to write to the fortran logical output */
/*     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that */
/*     this  logical  unit  must  be  attached  to  a  file or terminal */
/*     before calling this routine with a non-zero  value  for   IUNIT. */
/*     This routine does not check for the validity of a non-zero IUNIT */
/*     unit number. */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  MATVEC, MTTVEC, MSOLVE and the BLAS */
/* ***COMMON BLOCKS    SOLBLK */
/* ***END PROLOGUE  ISSCGN */

/* ***FIRST EXECUTABLE STATEMENT  ISSCGN */
    /* Parameter adjustments */
    --atdz;
    --dz;
    --atz;
    --atp;
    --p;
    --z__;
    --r__;
    --a;
    --x;
    --b;
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
	    (*msolve)(n, &b[1], &dz[1], nelt, &ia[1], &ja[1], &a[1], isym, &
		    rwork[1], &iwork[1]);
	    (*mttvec)(n, &dz[1], &atdz[1], nelt, &ia[1], &ja[1], &a[1], isym);
	    *bnrm = snrm2_(n, &atdz[1], &c__1);
	}
	*err = snrm2_(n, &atz[1], &c__1) / *bnrm;
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
	    io___46.ciunit = *iunit;
	    s_wsfe(&io___46);
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*itol), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	io___47.ciunit = *iunit;
	s_wsfe(&io___47);
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
/* ------------- LAST LINE OF ISSCGN FOLLOWS ---------------------------- */
} /* isscgn_ */

