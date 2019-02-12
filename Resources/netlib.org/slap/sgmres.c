/* sgmres.f -- translated by f2c (version 20100827).
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
static integer c__20 = 20;

/* DECK SGMRES */
/* Subroutine */ int sgmres_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, S_fp matvec, S_fp 
	msolve, integer *itol, real *tol, integer *itmax, integer *iter, real 
	*err, integer *ierr, integer *iunit, real *sb, real *sx, real *rgwk, 
	integer *lrgw, integer *igwk, integer *ligw, real *rwork, integer *
	iwork)
{
    /* System generated locals */
    integer i__1;
    real r__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, lq, lr, lv, lw, lz, ldl, kmp, nms, lxl;
    static real sum;
    static integer lzm1, lhes;
    static real bnrm;
    static integer jpre, maxl, lgmr, nmsl;
    static real rhol;
    extern doublereal snrm2_(integer *, real *, integer *);
    static integer iflag, jscal, nrmax;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    static integer nrsts;
    extern doublereal r1mach_(integer *);
    static integer maxlp1;
    extern /* Subroutine */ int spigmr_(integer *, real *, real *, real *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    S_fp, S_fp, integer *, real *, real *, real *, real *, integer *, 
	    real *, integer *, real *, real *, real *, integer *, real *, 
	    real *, real *, real *, integer *, real *, integer *, integer *, 
	    integer *, real *, integer *, integer *, integer *, real *);

/* ***BEGIN PROLOGUE  SGMRES */
/* ***DATE WRITTEN   871001   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SGMRES-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Generalized Minimum Residual */
/* ***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov */
/*           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/* ***PURPOSE  Preconditioned GMRES iterative sparse Ax=b solver. */
/*            This routine uses the generalized minimum residual */
/*            (GMRES) method with preconditioning to solve */
/*            non-symmetric linear systems of the form: A*x = b. */
/* ***DESCRIPTION */
/* *Usage: */
/*      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX */
/*      INTEGER   IERR, IUNIT, LRGW, LIGW, IGWK(LIGW) */
/*      INTEGER   IWORK(USER DEFINED) */
/*      REAL      B(N), X(N), A(NELT), TOL, ERR, SB(N), SX(N), */
/*      REAL      RGWK(LRGW), RWORK(USER DEFINED) */
/*      EXTERNAL  MATVEC, MSOLVE */

/*      CALL SGMRES(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, */
/*     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, SB, SX, */
/*     $     RGWK, LRGW, IGWK, LIGW, RWORK, IWORK) */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right-hand side vector. */
/* X      :INOUT    Real X(N). */
/*         On input X is your initial guess for the solution vector. */
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
/*         Y = A*X given A and X.  The name of the MATVEC routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MATVEC is: */
/*             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM ) */
/*         where N is the number of unknowns, Y is the product A*X */
/*         upon return, X is an input vector, and NELT is the number of */
/*         non-zeros in the SLAP IA, JA, A storage for the matrix A. */
/*         ISYM is a flag which, if non-zero, denotes that A is */
/*         symmetric and only the lower or upper triangle is stored. */
/* MSOLVE :EXT      External. */
/*         Name of the routine which solves a linear system Mz = r for */
/*         z given r with the preconditioning matrix M (M is supplied via */
/*         RWORK and IWORK arrays.  The name of the MSOLVE routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MSLOVE is: */
/*             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK) */
/*         Where N is the number of unknowns, R is the right-hand side */
/*         vector, and z is the solution upon return.  RWORK is a real */
/*         array that can be used to pass necessary preconditioning */
/*         information and/or workspace to MSOLVE.  IWORK is an integer */
/*         work array for the same purpose as RWORK. */
/* ITOL   :IN       Integer. */
/*         Flag to indicate the type of convergence criterion used. */
/*         ITOL=0  Means the  iteration stops when the test described */
/*                 below on  the  residual RL  is satisfied.  This is */
/*                 the  "Natural Stopping Criteria" for this routine. */
/*                 Other values  of   ITOL  cause  extra,   otherwise */
/*                 unnecessary, computation per iteration and     are */
/*                 therefore  much less  efficient.  See  ISSGMR (the */
/*                 stop test routine) for more information. */
/*         ITOL=1  Means   the  iteration stops   when the first test */
/*                 described below on  the residual RL  is satisfied, */
/*                 and there  is either right  or  no preconditioning */
/*                 being used. */
/*         ITOL=2  Implies     that   the  user    is   using    left */
/*                 preconditioning, and the second stopping criterion */
/*                 below is used. */
/*         ITOL=3  Means the  iteration stops   when  the  third test */
/*                 described below on Minv*Residual is satisfied, and */
/*                 there is either left  or no  preconditioning begin */
/*                 used. */
/*         ITOL=11 is    often  useful  for   checking  and comparing */
/*                 different routines.  For this case, the  user must */
/*                 supply  the  "exact" solution or  a  very accurate */
/*                 approximation (one with  an  error much less  than */
/*                 TOL) through a common block, */
/*                     COMMON /SOLBLK/ SOLN(1) */
/*                 if ITOL=11, iteration stops when the 2-norm of the */
/*                 difference between the iterative approximation and */
/*                 the user-supplied solution  divided by the  2-norm */
/*                 of the  user-supplied solution  is  less than TOL. */
/*                 Note that this requires  the  user to  set up  the */
/*                 "COMMON     /SOLBLK/ SOLN(LENGTH)"  in the calling */
/*                 routine.  The routine with this declaration should */
/*                 be loaded before the stop test so that the correct */
/*                 length is used by  the loader.  This procedure  is */
/*                 not standard Fortran and may not work correctly on */
/*                 your   system (although  it  has  worked  on every */
/*                 system the authors have tried).  If ITOL is not 11 */
/*                 then this common block is indeed standard Fortran. */
/* TOL    :INOUT    Real. */
/*         Convergence criterion, as described below.  If TOL is set */
/*         to zero on input, then a default value of 500*(the smallest */
/*         positive magnitude, machine epsilon) is used. */
/* ITMAX  :DUMMY    Integer. */
/*         Maximum number of iterations in most SLAP routines.  In */
/*         this routine this does not make sense.  The maximum number */
/*         of iterations here is given by ITMAX = MAXL*(NRMAX+1). */
/*         See IGWK for definitions of MAXL and NRMAX. */
/* ITER   :OUT      Integer. */
/*         Number of iterations required to reach convergence, or */
/*         ITMAX if convergence criterion could not be achieved in */
/*         ITMAX iterations. */
/* ERR    :OUT      Real. */
/*         Error estimate of error in final approximate solution, as */
/*         defined by ITOL.  Letting norm() denote the Euclidean */
/*         norm, ERR is defined as follows.. */

/*         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B), */
/*                               for right or no preconditioning, and */
/*                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/ */
/*                                norm(SB*(M-inverse)*B), */
/*                               for left preconditioning. */
/*         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B), */
/*                               since right or no preconditioning */
/*                               being used. */
/*         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/ */
/*                                norm(SB*(M-inverse)*B), */
/*                               since left preconditioning is being */
/*                               used. */
/*         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)| */
/*                               i=1,n */
/*         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN). */
/* IERR   :OUT      Integer. */
/*         Return error flag. */
/*               IERR = 0 => All went well. */
/*               IERR = 1 => Insufficient storage allocated for */
/*                           RGWK or IGWK. */
/*               IERR = 2 => Routine SGMREs failed to reduce the norm */
/*                           of the current residual on its last call, */
/*                           and so the iteration has stalled.  In */
/*                           this case, X equals the last computed */
/*                           approximation.  The user must either */
/*                           increase MAXL, or choose a different */
/*                           initial guess. */
/*               IERR =-1 => Insufficient length for RGWK array. */
/*                           IGWK(6) contains the required minimum */
/*                           length of the RGWK array. */
/*               IERR =-2 => Inconsistent ITOL and JPRE values. */
/*         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the */
/*         left-hand-side of the relevant stopping test defined */
/*         below associated with the residual for the current */
/*         approximation X(L). */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* SB     :IN       Real SB(N). */
/*         Array of length N containing scale factors for the right */
/*         hand side vector B.  If JSCAL.eq.0 (see below), SB need */
/*         not be supplied. */
/* SX     :IN       Real SX(N). */
/*         Array of length N containing scale factors for the solution */
/*         vector X.  If JSCAL.eq.0 (see below), SX need not be */
/*         supplied.  SB and SX can be the same array in the calling */
/*         program if desired. */
/* RGWK   :INOUT    Real RGWK(LRGW). */
/*         Real array of size at least 1 + N*(MAXL+6) + MAXL*(MAXL+3) */
/*         used for work space by SGMRES.  See below for definition of */
/*         MAXL. */
/*         On return, RGWK(1) = RHOL.  See IERR for definition of RHOL. */
/* LRGW   :IN       Integer. */
/*         Length of the real workspace, RGWK.  LRGW > 1 + N*(MAXL+6) */
/*         + MAXL*(MAXL+3). For the default values, RGWK has size at */
/*         least 131 + 16*N. */
/* IGWK   :INOUT    Integer IGWK(LIGW). */
/*         The following IGWK parameters should be set by the user */
/*         before calling this routine. */
/*         IGWK(1) = MAXL.  Maximum dimension of Krylov subspace in */
/*            which X - X0 is to be found (where, X0 is the initial */
/*            guess).  The default value of MAXL is 10. */
/*         IGWK(2) = KMP.  Maximum number of previous Krylov basis */
/*            vectors to which each new basis vector is made orthogonal. */
/*            The default value of KMP is MAXL. */
/*         IGWK(3) = JSCAL.  Flag indicating whether the scaling */
/*            arrays SB and SX are to be used. */
/*            JSCAL = 0 => SB and SX are not used and the algorithm */
/*               will perform as if all SB(I) = 1 and SX(I) = 1. */
/*            JSCAL = 1 =>  Only SX is used, and the algorithm */
/*               performs as if all SB(I) = 1. */
/*            JSCAL = 2 =>  Only SB is used, and the algorithm */
/*               performs as if all SX(I) = 1. */
/*            JSCAL = 3 =>  Both SB and SX are used. */
/*         IGWK(4) = JPRE.  Flag indicating whether preconditioning */
/*            is being used. */
/*            JPRE = 0  =>  There is no preconditioning. */
/*            JPRE > 0  =>  There is preconditioning on the right */
/*               only, and the solver will call routine MSOLVE. */
/*            JPRE < 0  =>  There is preconditioning on the left */
/*               only, and the solver will call routine MSOLVE. */
/*         IGWK(5) = NRMAX.  Maximum number of restarts of the */
/*            Krylov iteration.  The default value of NRMAX = 10. */
/*            if IWORK(5) = -1,  then no restarts are performed (in */
/*            this case, NRMAX is set to zero internally). */
/*         The following IWORK parameters are diagnostic information */
/*         made available to the user after this routine completes. */
/*         IGWK(6) = MLWK.  Required minimum length of RGWK array. */
/*         IGWK(7) = NMS.  The total number of calls to MSOLVE. */
/* LIGW   :IN       Integer. */
/*         Length of the integer workspace, IGWK.  LIGW >= 20. */

/* *Description: */
/*       SGMRES solves a linear system A*X = B rewritten in the form: */

/*        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B, */

/*       with right preconditioning, or */

/*        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B, */

/*       with left preconditioning, where A is an N-by-N real matrix, */
/*       X  and  B are N-vectors,   SB and SX   are  diagonal scaling */
/*       matrices,   and M is  a preconditioning    matrix.   It uses */
/*       preconditioned  Krylov   subpace  methods  based     on  the */
/*       generalized minimum residual  method (GMRES).   This routine */
/*       optionally performs  either  the  full     orthogonalization */
/*       version of the  GMRES  algorithm or an incomplete variant of */
/*       it.  Both versions use restarting of the linear iteration by */
/*       default, although the user can disable this feature. */

/*       The GMRES  algorithm generates a sequence  of approximations */
/*       X(L) to the  true solution of the above  linear system.  The */
/*       convergence criteria for stopping the  iteration is based on */
/*       the size  of the  scaled norm of  the residual  R(L)  =  B - */
/*       A*X(L).  The actual stopping test is either: */

/*               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B), */

/*       for right preconditioning, or */

/*               norm(SB*(M-inverse)*(B-A*X(L))) .le. */
/*                       TOL*norm(SB*(M-inverse)*B), */

/*       for left preconditioning, where norm() denotes the euclidean */
/*       norm, and TOL is  a positive scalar less  than one  input by */
/*       the user.  If TOL equals zero  when SGMRES is called, then a */
/*       default  value  of 500*(the   smallest  positive  magnitude, */
/*       machine epsilon) is used.  If the  scaling arrays SB  and SX */
/*       are used, then  ideally they  should be chosen  so  that the */
/*       vectors SX*X(or SX*M*X) and  SB*B have all their  components */
/*       approximately equal  to  one in  magnitude.  If one wants to */
/*       use the same scaling in X  and B, then  SB and SX can be the */
/*       same array in the calling program. */

/*       The following is a list of the other routines and their */
/*       functions used by SGMRES: */
/*       SPIGMR  Contains the main iteration loop for GMRES. */
/*       SORTH   Orthogonalizes a new vector against older basis vects. */
/*       SHEQR   Computes a QR decomposition of a Hessenberg matrix. */
/*       SHELS   Solves a Hessenberg least-squares system, using QR */
/*               factors. */
/*       SRLCAL  Computes the scaled residual RL. */
/*       SXLCAL  Computes the solution XL. */
/*       ISSGMR  User-replaceable stopping routine. */

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
/* ***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh, */
/*                 "Reduced Storage Matrix Methods In Stiff ODE */
/*                 Systems," LLNL report UCRL-95088, Rev. 1, */
/*                 June 1987. */
/* ***ROUTINES CALLED  SPIGMR, SORTHO, SHEQR, SHELS, SRCAL, SXLCAL, */
/*                    ISSGMR, SNRM2, SDOT, SAXPY, SSCAL, ISAMAX, R1MACH. */
/* ***END PROLOGUE  SGMRES */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */

/* ***FIRST EXECUTABLE STATEMENT  SGMRES */
    /* Parameter adjustments */
    --sx;
    --sb;
    --x;
    --b;
    --a;
    --ja;
    --ia;
    --rgwk;
    --igwk;
    --rwork;
    --iwork;

    /* Function Body */
    *ierr = 0;
/*   ------------------------------------------------------------------ */
/*         Load method parameters with user values or defaults. */
/*   ------------------------------------------------------------------ */
    maxl = igwk[1];
    if (maxl == 0) {
	maxl = 10;
    }
    if (maxl > *n) {
	maxl = *n;
    }
    kmp = igwk[2];
    if (kmp == 0) {
	kmp = maxl;
    }
    if (kmp > maxl) {
	kmp = maxl;
    }
    jscal = igwk[3];
    jpre = igwk[4];
/*         Check for consistent values of ITOL and JPRE. */
    if (*itol == 1 && jpre < 0) {
	goto L650;
    }
    if (*itol == 2 && jpre >= 0) {
	goto L650;
    }
    nrmax = igwk[5];
    if (nrmax == 0) {
	nrmax = 10;
    }
/*         If NRMAX .eq. -1, then set NRMAX = 0 to turn off restarting. */
    if (nrmax == -1) {
	nrmax = 0;
    }
/*         If input value of TOL is zero, set it to its default value. */
    if (*tol == 0.f) {
	*tol = r1mach_(&c__3) * 500.f;
    }

/*         Initialize counters. */
    *iter = 0;
    nms = 0;
    nrsts = 0;
/*   ------------------------------------------------------------------ */
/*         Form work array segment pointers. */
/*   ------------------------------------------------------------------ */
    maxlp1 = maxl + 1;
    lv = 1;
    lr = lv + *n * maxlp1;
    lhes = lr + *n + 1;
    lq = lhes + maxl * maxlp1;
    ldl = lq + (maxl << 1);
    lw = ldl + *n;
    lxl = lw + *n;
    lz = lxl + *n;

/*         Load igwk(6) with required minimum length of the rgwk array. */
    igwk[6] = lz + *n - 1;
    if (lz + *n - 1 > *lrgw) {
	goto L640;
    }
/*   ------------------------------------------------------------------ */
/*         Calculate scaled-preconditioned norm of RHS vector b. */
/*   ------------------------------------------------------------------ */
    if (jpre < 0) {
	(*msolve)(n, &b[1], &rgwk[lr], nelt, &ia[1], &ja[1], &a[1], isym, &
		rwork[1], &iwork[1]);
	++nms;
    } else {
	scopy_(n, &b[1], &c__1, &rgwk[lr], &c__1);
    }
    if (jscal == 2 || jscal == 3) {
	sum = 0.f;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    r__1 = rgwk[lr - 1 + i__] * sb[i__];
	    sum += r__1 * r__1;
/* L10: */
	}
	bnrm = sqrt(sum);
    } else {
	bnrm = snrm2_(n, &rgwk[lr], &c__1);
    }
/*   ------------------------------------------------------------------ */
/*         Calculate initial residual. */
/*   ------------------------------------------------------------------ */
    (*matvec)(n, &x[1], &rgwk[lr], nelt, &ia[1], &ja[1], &a[1], isym);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rgwk[lr - 1 + i__] = b[i__] - rgwk[lr - 1 + i__];
/* L50: */
    }
/*   ------------------------------------------------------------------ */
/*         If performing restarting, then load the residual into the */
/*         correct location in the Rgwk array. */
/*   ------------------------------------------------------------------ */
L100:
    if (nrsts > nrmax) {
	goto L610;
    }
    if (nrsts > 0) {
/*         Copy the curr residual to different loc in the Rgwk array. */
	scopy_(n, &rgwk[ldl], &c__1, &rgwk[lr], &c__1);
    }
/*   ------------------------------------------------------------------ */
/*         Use the SPIGMR algorithm to solve the linear system A*Z = R. */
/*   ------------------------------------------------------------------ */
    spigmr_(n, &rgwk[lr], &sb[1], &sx[1], &jscal, &maxl, &maxlp1, &kmp, &
	    nrsts, &jpre, (S_fp)matvec, (S_fp)msolve, &nmsl, &rgwk[lz], &rgwk[
	    lv], &rgwk[lhes], &rgwk[lq], &lgmr, &rwork[1], &iwork[1], &rgwk[
	    lw], &rgwk[ldl], &rhol, &nrmax, &b[1], &bnrm, &x[1], &rgwk[lxl], 
	    itol, tol, nelt, &ia[1], &ja[1], &a[1], isym, iunit, &iflag, err);
    *iter += lgmr;
    nms += nmsl;

/*         Increment X by the current approximate solution Z of A*Z = R. */

    lzm1 = lz - 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] += rgwk[lzm1 + i__];
/* L110: */
    }
    if (iflag == 0) {
	goto L600;
    }
    if (iflag == 1) {
	++nrsts;
	goto L100;
    }
    if (iflag == 2) {
	goto L620;
    }
/*   ------------------------------------------------------------------ */
/*         All returns are made through this section. */
/*   ------------------------------------------------------------------ */
/*         The iteration has converged. */

L600:
    igwk[7] = nms;
    rgwk[1] = rhol;
    *ierr = 0;
    return 0;

/*         Max number((NRMAX+1)*MAXL) of linear iterations performed. */
L610:
    igwk[7] = nms;
    rgwk[1] = rhol;
    *ierr = 1;
    return 0;

/*         GMRES failed to reduce last residual in MAXL iterations. */
/*         The iteration has stalled. */
L620:
    igwk[7] = nms;
    rgwk[1] = rhol;
    *ierr = 2;
    return 0;
/*         Error return.  Insufficient length for Rgwk array. */
L640:
    *err = *tol;
    *ierr = -1;
    return 0;
/*         Error return.  Inconsistent ITOL and JPRE values. */
L650:
    *err = *tol;
    *ierr = -2;
    return 0;
/* ------------- LAST LINE OF SGMRES FOLLOWS ---------------------------- */
} /* sgmres_ */

/* DECK SSDGMR */
/* Subroutine */ int ssdgmr_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, integer *nsave, 
	integer *itol, real *tol, integer *itmax, integer *iter, real *err, 
	integer *ierr, integer *iunit, real *rwork, integer *lenw, integer *
	iwork, integer *leniw)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int ss2y_(integer *, integer *, integer *, 
	    integer *, real *, integer *), ssdi_();
    static integer locw;
    extern /* Subroutine */ int ssds_(integer *, integer *, integer *, 
	    integer *, real *, integer *, real *), ssmv_();
    static integer lociw;
    extern /* Subroutine */ int schkw_(char *, integer *, integer *, integer *
	    , integer *, integer *, integer *, real *, ftnlen);
    static integer locdin, locigw, locrgw;
    extern /* Subroutine */ int sgmres_(integer *, real *, real *, integer *, 
	    integer *, integer *, real *, integer *, S_fp, S_fp, integer *, 
	    real *, integer *, integer *, real *, integer *, integer *, real *
	    , real *, real *, integer *, integer *, integer *, real *, 
	    integer *);
    static integer myitol;

/* ***BEGIN PROLOGUE  SSDGMR */
/* ***DATE WRITTEN   880615   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSDGMR-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Generalized Minimum Residual */
/* ***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov */
/*           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/* ***PURPOSE  Diagonally scaled GMRES iterative sparse Ax=b solver. */
/*            This routine uses the generalized minimum residual */
/*            (GMRES) method with diagonal scaling to solve possibly */
/*            non-symmetric linear systems of the form: A*x = b. */
/* ***DESCRIPTION */
/* *Usage: */
/*      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE */
/*      INTEGER   ITOL, ITMAX, IERR, IUNIT, LENW, IWORK(LENIW), LENIW */
/*      REAL      B(N), X(N), A(NELT), TOL, ERR */
/*      REAL      RWORK(LENW) */
/*      EXTERNAL  MATVEC, MSOLVE */

/*      CALL SSDGMR(N, B, X, NELT, IA, JA, A, ISYM, NSAVE, */
/*     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, */
/*     $     RWORK, LENW, IWORK, LENIW) */

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
/*         These arrays should hold the matrix A in either the SLAP */
/*         Triad format or the SLAP Column format.  See "Description", */
/*         below.  If the SLAP Triad format is chosen it is changed */
/*         internally to the SLAP Column format. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* NSAVE  :IN       Integer. */
/*         Number of direction vectors to save and orthogonalize against. */
/*         Must be greater than 1. */
/* ITOL   :IN       Integer. */
/*         Flag to indicate the type of convergence criterion used. */
/*         ITOL=0  Means the  iteration stops when the test described */
/*                 below on  the  residual RL  is satisfied.  This is */
/*                 the  "Natural Stopping Criteria" for this routine. */
/*                 Other values  of   ITOL  cause  extra,   otherwise */
/*                 unnecessary, computation per iteration and     are */
/*                 therefore  much less  efficient.  See  ISSGMR (the */
/*                 stop test routine) for more information. */
/*         ITOL=1  Means   the  iteration stops   when the first test */
/*                 described below on  the residual RL  is satisfied, */
/*                 and there  is either right  or  no preconditioning */
/*                 being used. */
/*         ITOL=2  Implies     that   the  user    is   using    left */
/*                 preconditioning, and the second stopping criterion */
/*                 below is used. */
/*         ITOL=3  Means the  iteration stops   when  the  third test */
/*                 described below on Minv*Residual is satisfied, and */
/*                 there is either left  or no  preconditioning begin */
/*                 used. */
/*         ITOL=11 is    often  useful  for   checking  and comparing */
/*                 different routines.  For this case, the  user must */
/*                 supply  the  "exact" solution or  a  very accurate */
/*                 approximation (one with  an  error much less  than */
/*                 TOL) through a common block, */
/*                     COMMON /SOLBLK/ SOLN(1) */
/*                 if ITOL=11, iteration stops when the 2-norm of the */
/*                 difference between the iterative approximation and */
/*                 the user-supplied solution  divided by the  2-norm */
/*                 of the  user-supplied solution  is  less than TOL. */
/*                 Note that this requires  the  user to  set up  the */
/*                 "COMMON     /SOLBLK/ SOLN(LENGTH)"  in the calling */
/*                 routine.  The routine with this declaration should */
/*                 be loaded before the stop test so that the correct */
/*                 length is used by  the loader.  This procedure  is */
/*                 not standard Fortran and may not work correctly on */
/*                 your   system (although  it  has  worked  on every */
/*                 system the authors have tried).  If ITOL is not 11 */
/*                 then this common block is indeed standard Fortran. */
/* TOL    :INOUT    Real. */
/*         Convergence criterion, as described below.  If TOL is set */
/*         to zero on input, then a default value of 500*(the smallest */
/*         positive magnitude, machine epsilon) is used. */
/* ITMAX  :IN       Integer. */
/*         Maximum number of iterations.  This routine uses the default */
/*         of NRMAX = ITMAX/NSAVE to determine the when each restart */
/*         oshould ccur.  See the description of NRMAX and MAXL in */
/*         SGMRES for a full and frightfully interesting discussion of */
/*         this topic. */
/* ITER   :OUT      Integer. */
/*         Number of iterations required to reach convergence, or */
/*         ITMAX+1 if convergence criterion could not be achieved in */
/*         ITMAX iterations. */
/* ERR    :OUT      Real. */
/*         Error estimate of error in final approximate solution, as */
/*         defined by ITOL.  Letting norm() denote the Euclidean */
/*         norm, ERR is defined as follows... */
/*         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B), */
/*                               for right or no preconditioning, and */
/*                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/ */
/*                                norm(SB*(M-inverse)*B), */
/*                               for left preconditioning. */
/*         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B), */
/*                               since right or no preconditioning */
/*                               being used. */
/*         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/ */
/*                                norm(SB*(M-inverse)*B), */
/*                               since left preconditioning is being */
/*                               used. */
/*         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)| */
/*                               i=1,n */
/*         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN). */
/* IERR   :OUT      Integer. */
/*         Return error flag. */
/*               IERR = 0 => All went well. */
/*               IERR = 1 => Insufficient storage allocated for */
/*                           RGWK or IGWK. */
/*               IERR = 2 => Routine SPIGMR failed to reduce the norm */
/*                           of the current residual on its last call, */
/*                           and so the iteration has stalled.  In */
/*                           this case, X equals the last computed */
/*                           approximation.  The user must either */
/*                           increase MAXL, or choose a different */
/*                           initial guess. */
/*               IERR =-1 => Insufficient length for RGWK array. */
/*                           IGWK(6) contains the required minimum */
/*                           length of the RGWK array. */
/*               IERR =-2 => Inconsistent ITOL and JPRE values. */
/*         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the */
/*         left-hand-side of the relevant stopping test defined */
/*         below associated with the residual for the current */
/*         approximation X(L). */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* RWORK  :WORK    Real RWORK(LENW). */
/*         Real array of size LENW. */
/* LENW   :IN       Integer. */
/*         Length of the real workspace, RWORK.  LENW >= 1 + N*(NSAVE+7) */
/*         + NSAVE*(NSAVE+3). For the recommended values of NSAVE (10), */
/*         RWORK has size at least 131 + 17*N. */
/* IWORK  :INOUT    Integer IWORK(USER DEFINED >= 30). */
/*         Used to hold pointers into the RWORK array. */
/*         Upon return the following locations of IWORK hold information */
/*         which may be of use to the user: */
/*         IWORK(9)  Amount of Integer workspace actually used. */
/*         IWORK(10) Amount of Real    workspace actually used. */
/* LENIW  :IN       Integer. */
/*         Length of the integer workspace IWORK.  LENIW >= 30. */

/* *Description: */
/*       SSDGMR solves a linear system A*X = B rewritten in the form: */

/*        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B, */

/*       with right preconditioning, or */

/*        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B, */

/*       with left preconditioning, where a is an n-by-n real matrix, */
/*       X and  B  are N-vectors,  SB and  SX  are  diagonal  scaling */
/*       matrices, and  M   is   the  diagonal  of   A.     It   uses */
/*       preconditioned   Krylov  subpace   methods  based    on  the */
/*       generalized  minimum residual method (GMRES).   This routine */
/*       is  a  driver routine  which   assumes a  SLAP matrix   data */
/*       structure  and   sets  up the  necessary information   to do */
/*       diagonal preconditioning and  calls  the main GMRES  routine */
/*       SGMRES   for  the  solution  of the   linear system.  SGMRES */
/*       optionally   performs   either the   full  orthogonalization */
/*       version of the GMRES algorithm or an  incomplete  variant of */
/*       it.  Both versions use restarting of the linear iteration by */
/*       default, although the user can disable this feature. */

/*       The GMRES  algorithm generates a sequence  of approximations */
/*       X(L) to the  true solution of the above  linear system.  The */
/*       convergence criteria for stopping the  iteration is based on */
/*       the size  of the  scaled norm of  the residual  R(L)  =  B - */
/*       A*X(L).  The actual stopping test is either: */

/*               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B), */

/*       for right preconditioning, or */

/*               norm(SB*(M-inverse)*(B-A*X(L))) .le. */
/*                       TOL*norm(SB*(M-inverse)*B), */

/*       for left preconditioning, where norm() denotes the euclidean */
/*       norm, and TOL is  a positive scalar less  than one  input by */
/*       the user.  If TOL equals zero  when SSDGMR is called, then a */
/*       default  value  of 500*(the   smallest  positive  magnitude, */
/*       machine epsilon) is used.  If the  scaling arrays SB  and SX */
/*       are used, then  ideally they  should be chosen  so  that the */
/*       vectors SX*X(or SX*M*X) and  SB*B have all their  components */
/*       approximately equal  to  one in  magnitude.  If one wants to */
/*       use the same scaling in X  and B, then  SB and SX can be the */
/*       same array in the calling program. */

/*       The following is a list of the other routines and their */
/*       functions used by GMRES: */
/*       SGMRES  Contains the matrix structure independent driver */
/*               routine for GMRES. */
/*       SPIGMR  Contains the main iteration loop for GMRES. */
/*       SORTH   Orthogonalizes a new vector against older basis vects. */
/*       SHEQR   Computes a QR decomposition of a Hessenberg matrix. */
/*       SHELS   Solves a Hessenberg least-squares system, using QR */
/*               factors. */
/*       RLCALC  Computes the scaled residual RL. */
/*       XLCALC  Computes the solution XL. */
/*       ISSGMR  User-replaceable stopping routine. */

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
/* ***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh, */
/*                 "Reduced Storage Matrix Methods In Stiff ODE */
/*                 Systems," LLNL report UCRL-95088, Rev. 1, */
/*                 June 1987. */
/* ***ROUTINES CALLED  SS2Y, SCHKW, SSDS, SGMRES */
/* ***END PROLOGUE  SSDGMR */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */

/*         Change the SLAP input matrix IA, JA, A to SLAP-Column format. */
/* ***FIRST EXECUTABLE STATEMENT  SSDGMR */
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
    *err = 0.f;
    if (*nsave <= 1) {
	*ierr = 3;
	return 0;
    }
    ss2y_(n, nelt, &ia[1], &ja[1], &a[1], isym);

/*         Set up the workspace.  We assume MAXL=KMP=NSAVE. */
/*         Compute the inverse of the diagonal of the matrix. */
    locigw = 11;
    lociw = locigw + 20;

    locdin = 1;
    locrgw = locdin + *n;
    locw = locrgw + 1 + *n * (*nsave + 6) + *nsave * (*nsave + 3);

    iwork[4] = locdin;
    iwork[9] = lociw;
    iwork[10] = locw;

/*         Check the workspace allocations. */
    schkw_("SSDGMR", &lociw, leniw, &locw, lenw, ierr, iter, err, (ftnlen)6);
    if (*ierr != 0) {
	return 0;
    }

    ssds_(n, nelt, &ia[1], &ja[1], &a[1], isym, &rwork[locdin]);

/*         Perform the Diagonaly Scaled Generalized Minimum */
/*         Residual iteration algorithm.  The following SGMRES */
/*         defaults are used MAXL = KMP = NSAVE, JSCAL = 0, */
/*         JPRE = -1, NRMAX = ITMAX/NSAVE */
    iwork[locigw] = *nsave;
    iwork[locigw + 1] = *nsave;
    iwork[locigw + 2] = 0;
    iwork[locigw + 3] = -1;
    iwork[locigw + 4] = *itmax / *nsave;
    myitol = 0;

    i__1 = *lenw - locrgw;
    sgmres_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)ssmv_, (
	    S_fp)ssdi_, &myitol, tol, itmax, iter, err, ierr, iunit, &rwork[1]
	    , &rwork[1], &rwork[locrgw], &i__1, &iwork[locigw], &c__20, &
	    rwork[1], &iwork[1]);

    if (*iter > *itmax) {
	*ierr = 2;
    }
    return 0;
/* ------------- LAST LINE OF SSDGMR FOLLOWS ---------------------------- */
} /* ssdgmr_ */

/* DECK SSLUGM */
/* Subroutine */ int sslugm_(integer *n, real *b, real *x, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, integer *nsave, 
	integer *itol, real *tol, integer *itmax, integer *iter, real *err, 
	integer *ierr, integer *iunit, real *rwork, integer *lenw, integer *
	iwork, integer *leniw)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer j, nl, nu;
    extern /* Subroutine */ int ss2y_(integer *, integer *, integer *, 
	    integer *, real *, integer *);
    static integer jbgn, jend, icol, locl, locu, locw;
    extern /* Subroutine */ int ssmv_();
    static integer locnc, locil, locjl, lociu, locju, locnr, lociw;
    extern /* Subroutine */ int schkw_(char *, integer *, integer *, integer *
	    , integer *, integer *, integer *, real *, ftnlen), sslui_();
    static integer locdin, locigw, locrgw;
    extern /* Subroutine */ int sgmres_(integer *, real *, real *, integer *, 
	    integer *, integer *, real *, integer *, S_fp, S_fp, integer *, 
	    real *, integer *, integer *, real *, integer *, integer *, real *
	    , real *, real *, integer *, integer *, integer *, real *, 
	    integer *);
    static integer myitol;
    extern /* Subroutine */ int ssilus_(integer *, integer *, integer *, 
	    integer *, real *, integer *, integer *, integer *, integer *, 
	    real *, real *, integer *, integer *, integer *, real *, integer *
	    , integer *);

/* ***BEGIN PROLOGUE  SSLUGM */
/* ***DATE WRITTEN   880615   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SSLUGM-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Generalized Minimum Residual */
/* ***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov */
/*           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/* ***PURPOSE  Incomplete LU GMRES iterative sparse Ax=b solver. */
/*            This routine uses the generalized minimum residual */
/*            (GMRES) method with incomplete LU factorization for */
/*            preconditioning to solve possibly non-symmetric linear */
/*            systems of the form: Ax = b. */
/* ***DESCRIPTION */
/* *Usage: */
/*      INTEGER   N, NELT, IA(NELT), JA(NELT), ISYM, NSAVE */
/*      INTEGER   ITOL, ITMAX, IERR, IUNIT, LENW, IWORK(LENIW), LENIW */
/*      REAL      B(N), X(N), A(NELT), TOL, ERR, SB(N), SX(N), */
/*      REAL      RWORK(LENW) */
/*      EXTERNAL  MATVEC, MSOLVE */

/*      CALL SSLUGM(N, B, X, NELT, IA, JA, A, ISYM, NSAVE, */
/*     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, */
/*     $     RWORK, LENW, IWORK, LENIW) */

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
/*         These arrays should hold the matrix A in either the SLAP */
/*         Triad format or the SLAP Column format.  See "Description", */
/*         below.  If the SLAP Triad format is chosen it is changed */
/*         internally to the SLAP Column format. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* NSAVE  :IN       Integer. */
/*         Number of direction vectors to save and orthogonalize against. */
/*         Must be greater than 1. */
/* ITOL   :IN       Integer. */
/*         Flag to indicate the type of convergence criterion used. */
/*         ITOL=0  Means the  iteration stops when the test described */
/*                 below on  the  residual RL  is satisfied.  This is */
/*                 the  "Natural Stopping Criteria" for this routine. */
/*                 Other values  of   ITOL  cause  extra,   otherwise */
/*                 unnecessary, computation per iteration and     are */
/*                 therefore  much less  efficient.  See  ISSGMR (the */
/*                 stop test routine) for more information. */
/*         ITOL=1  Means   the  iteration stops   when the first test */
/*                 described below on  the residual RL  is satisfied, */
/*                 and there  is either right  or  no preconditioning */
/*                 being used. */
/*         ITOL=2  Implies     that   the  user    is   using    left */
/*                 preconditioning, and the second stopping criterion */
/*                 below is used. */
/*         ITOL=3  Means the  iteration stops   when  the  third test */
/*                 described below on Minv*Residual is satisfied, and */
/*                 there is either left  or no  preconditioning begin */
/*                 used. */
/*         ITOL=11 is    often  useful  for   checking  and comparing */
/*                 different routines.  For this case, the  user must */
/*                 supply  the  "exact" solution or  a  very accurate */
/*                 approximation (one with  an  error much less  than */
/*                 TOL) through a common block, */
/*                     COMMON /SOLBLK/ SOLN(1) */
/*                 if ITOL=11, iteration stops when the 2-norm of the */
/*                 difference between the iterative approximation and */
/*                 the user-supplied solution  divided by the  2-norm */
/*                 of the  user-supplied solution  is  less than TOL. */
/*                 Note that this requires  the  user to  set up  the */
/*                 "COMMON     /SOLBLK/ SOLN(LENGTH)"  in the calling */
/*                 routine.  The routine with this declaration should */
/*                 be loaded before the stop test so that the correct */
/*                 length is used by  the loader.  This procedure  is */
/*                 not standard Fortran and may not work correctly on */
/*                 your   system (although  it  has  worked  on every */
/*                 system the authors have tried).  If ITOL is not 11 */
/*                 then this common block is indeed standard Fortran. */
/* TOL    :INOUT    Real. */
/*         Convergence criterion, as described below.  If TOL is set */
/*         to zero on input, then a default value of 500*(the smallest */
/*         positive magnitude, machine epsilon) is used. */
/* ITMAX  :IN       Integer. */
/*         Maximum number of iterations.  This routine uses the default */
/*         of NRMAX = ITMAX/NSAVE to determine the when each restart */
/*         should occur.  See the description of NRMAX and MAXL in */
/*         SGMRES for a full and frightfully interesting discussion of */
/*         this topic. */
/* ITER   :OUT      Integer. */
/*         Number of iterations required to reach convergence, or */
/*         ITMAX+1 if convergence criterion could not be achieved in */
/*         ITMAX iterations. */
/* ERR    :OUT      Real. */
/*         Error estimate of error in final approximate solution, as */
/*         defined by ITOL.  Letting norm() denote the Euclidean */
/*         norm, ERR is defined as follows... */
/*         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B), */
/*                               for right or no preconditioning, and */
/*                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/ */
/*                                norm(SB*(M-inverse)*B), */
/*                               for left preconditioning. */
/*         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B), */
/*                               since right or no preconditioning */
/*                               being used. */
/*         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/ */
/*                                norm(SB*(M-inverse)*B), */
/*                               since left preconditioning is being */
/*                               used. */
/*         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)| */
/*                               i=1,n */
/*         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN). */
/* IERR   :OUT      Integer. */
/*         Return error flag. */
/*               IERR = 0 => All went well. */
/*               IERR = 1 => Insufficient storage allocated for */
/*                           RGWK or IGWK. */
/*               IERR = 2 => Routine SPIGMR failed to reduce the norm */
/*                           of the current residual on its last call, */
/*                           and so the iteration has stalled.  In */
/*                           this case, X equals the last computed */
/*                           approximation.  The user must either */
/*                           increase MAXL, or choose a different */
/*                           initial guess. */
/*               IERR =-1 => Insufficient length for RGWK array. */
/*                           IGWK(6) contains the required minimum */
/*                           length of the RGWK array. */
/*               IERR =-2 => Inconsistent ITOL and JPRE values. */
/*         For IERR <= 2, RGWK(1) = RHOL, which is the norm on the */
/*         left-hand-side of the relevant stopping test defined */
/*         below associated with the residual for the current */
/*         approximation X(L). */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* RWORK  :WORK    Real RWORK(LENW). */
/*         Real array of size LENW. */
/* LENW   :IN       Integer. */
/*         Length of the real workspace, RWORK. LENW >= 1 + N*(NSAVE+7) */
/*         +  NSAVE*(NSAVE+3)+NEL+NU. For the recommended values,  RWORK */
/*         has size at least 131 + 17*N + NEL + NU.  Where  NEL is  the */
/*         number of non- zeros  in  the  lower triangle of  the matrix */
/*         (including the diagonal).  NU is the  number  of nonzeros in */
/*         the upper triangle of the matrix (including the diagonal). */
/* IWORK  :INOUT    Integer IWORK(LENIW). */
/*         Used to hold pointers into the RWORK array. */
/*         Upon return the following locations of IWORK hold information */
/*         which may be of use to the user: */
/*         IWORK(9)  Amount of Integer workspace actually used. */
/*         IWORK(10) Amount of Real    workspace actually used. */
/* LENIW  :IN       Integer. */
/*         Length of the integer workspace, IWORK. */
/*         LENIW >= NEL+NU+4*N+32. */

/* *Description: */
/*       SSLUGM solves a linear system A*X = B rewritten in the form: */

/*        (SB*A*(M-inverse)*(SX-inverse))*(SX*M*X) = SB*B, */

/*       with right preconditioning, or */

/*        (SB*(M-inverse)*A*(SX-inverse))*(SX*X) = SB*(M-inverse)*B, */

/*       with left preconditioning, where a is an n-by-n real matrix, */
/*       X and  B are  N-vectors,  SB and  SX  are   diagonal scaling */
/*       matrices, and M is the Incomplete LU factorization of A.  It */
/*       uses preconditioned  Krylov subpace   methods  based on  the */
/*       generalized minimum residual  method (GMRES).   This routine */
/*       is a  driver  routine  which  assumes a SLAP   matrix   data */
/*       structure   and  sets  up  the  necessary  information to do */
/*       diagonal  preconditioning  and calls the main GMRES  routine */
/*       SGMRES for the   solution   of the linear   system.   SGMRES */
/*       optionally   performs  either  the full    orthogonalization */
/*       version of the  GMRES algorithm or  an incomplete variant of */
/*       it.  Both versions use restarting of the linear iteration by */
/*       default, although the user can disable this feature. */

/*       The GMRES  algorithm generates a sequence  of approximations */
/*       X(L) to the  true solution of the above  linear system.  The */
/*       convergence criteria for stopping the  iteration is based on */
/*       the size  of the  scaled norm of  the residual  R(L)  =  B - */
/*       A*X(L).  The actual stopping test is either: */

/*               norm(SB*(B-A*X(L))) .le. TOL*norm(SB*B), */

/*       for right preconditioning, or */

/*               norm(SB*(M-inverse)*(B-A*X(L))) .le. */
/*                       TOL*norm(SB*(M-inverse)*B), */

/*       for left preconditioning, where norm() denotes the euclidean */
/*       norm, and TOL is  a positive scalar less  than one  input by */
/*       the user.  If TOL equals zero  when SSLUGM is called, then a */
/*       default  value  of 500*(the   smallest  positive  magnitude, */
/*       machine epsilon) is used.  If the  scaling arrays SB  and SX */
/*       are used, then  ideally they  should be chosen  so  that the */
/*       vectors SX*X(or SX*M*X) and  SB*B have all their  components */
/*       approximately equal  to  one in  magnitude.  If one wants to */
/*       use the same scaling in X  and B, then  SB and SX can be the */
/*       same array in the calling program. */

/*       The following is a list of the other routines and their */
/*       functions used by GMRES: */
/*       SGMRES  Contains the matrix structure independent driver */
/*               routine for GMRES. */
/*       SPIGMR  Contains the main iteration loop for GMRES. */
/*       SORTH   Orthogonalizes a new vector against older basis vects. */
/*       SHEQR   Computes a QR decomposition of a Hessenberg matrix. */
/*       SHELS   Solves a Hessenberg least-squares system, using QR */
/*               factors. */
/*       RLCALC  Computes the scaled residual RL. */
/*       XLCALC  Computes the solution XL. */
/*       ISSGMR  User-replaceable stopping routine. */

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
/* ***REFERENCES  1. Peter N. Brown and A. C. Hindmarsh, */
/*                 "Reduced Storage Matrix Methods In Stiff ODE */
/*                 Systems," LLNL report UCRL-95088, Rev. 1, */
/*                 June 1987. */
/* ***ROUTINES CALLED  SS2Y, SCHKW, SSILUS, SGMRES, SSMV, SSLUI */
/* ***END PROLOGUE  SSLUGM */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */

/*         Change the SLAP input matrix IA, JA, A to SLAP-Column format. */
/* ***FIRST EXECUTABLE STATEMENT  SSLUGM */
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
    *err = 0.f;
    if (*nsave <= 1) {
	*ierr = 3;
	return 0;
    }
    ss2y_(n, nelt, &ia[1], &ja[1], &a[1], isym);

/*         Count number of Non-Zero elements preconditioner ILU matrix. */
/*         Then set up the work arrays.  We assume MAXL=KMP=NSAVE. */
    nl = 0;
    nu = 0;
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
/*         Don't count diagonal. */
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

    locigw = 11;
    locil = locigw + 20;
    locjl = locil + *n + 1;
    lociu = locjl + nl;
    locju = lociu + nu;
    locnr = locju + *n + 1;
    locnc = locnr + *n;
    lociw = locnc + *n;

    locl = 1;
    locdin = locl + nl;
    locu = locdin + *n;
    locrgw = locu + nu;
    locw = locrgw + 1 + *n * (*nsave + 6) + *nsave * (*nsave + 3);

/*         Check the workspace allocations. */
    schkw_("SSLUGM", &lociw, leniw, &locw, lenw, ierr, iter, err, (ftnlen)6);
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

/*         Perform the Incomplet LU Preconditioned Generalized Minimum */
/*         Residual iteration algorithm.  The following SGMRES */
/*         defaults are used MAXL = KMP = NSAVE, JSCAL = 0, */
/*         JPRE = -1, NRMAX = ITMAX/NSAVE */
    iwork[locigw] = *nsave;
    iwork[locigw + 1] = *nsave;
    iwork[locigw + 2] = 0;
    iwork[locigw + 3] = -1;
    iwork[locigw + 4] = *itmax / *nsave;
    myitol = 0;

    i__1 = *lenw - locrgw;
    sgmres_(n, &b[1], &x[1], nelt, &ia[1], &ja[1], &a[1], isym, (S_fp)ssmv_, (
	    S_fp)sslui_, &myitol, tol, itmax, iter, err, ierr, iunit, &rwork[
	    1], &rwork[1], &rwork[locrgw], &i__1, &iwork[locigw], &c__20, &
	    rwork[1], &iwork[1]);

    if (*iter > *itmax) {
	*ierr = 2;
    }
    return 0;
/* ------------- LAST LINE OF SSLUGM FOLLOWS ---------------------------- */
} /* sslugm_ */

/* DECK SHELS */
/* Subroutine */ int shels_(real *a, integer *lda, integer *n, real *q, real *
	b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static real c__;
    static integer k;
    static real s, t, t1, t2;
    static integer kb, iq, kp1;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);

/* ***BEGIN PROLOGUE  SHEQR */
/* ***DATE WRITTEN   871001   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SHEQR-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Generalized Minimum Residual */
/* ***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov */
/*           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/* ***PURPOSE  Internal routine for SGMRES. */
/* ***DESCRIPTION */
/*        This routine is extraced from the LINPACK routine SGESL with */
/*        changes  due to  the fact  that  A is an  upper   Hessenberg */
/*        matrix. */

/*        SHELS solves the least squares problem: */

/*                   MIN(B-A*X,B-A*X) */

/*        using the factors computed by SHEQR. */

/* *Usage: */
/*      INTEGER LDA, N */
/*      REAL A(LDA,1), B(1), Q(1) */

/*      CALL SHELS(A, LDA, N, Q, B) */

/* *Arguments: */
/* A       :IN       Real A(LDA,N) */
/*          The output from SHEQR which contains the upper */
/*          triangular factor R in the QR decomposition of A. */
/* LDA     :IN       Integer */
/*          The leading dimension of the array A. */
/* N       :IN       Integer */
/*          A is originally an (N+1) by N matrix. */
/* Q       :IN       Real Q(2*N) */
/*          The coefficients of the N givens rotations */
/*          used in the QR factorization of A. */
/* B       :INOUT    Real B(N+1) */
/*          On input, B is the right hand side vector. */
/*          On output, B is the solution vector X. */
/* *See Also: */
/*         SGMRES */

/* ***ROUTINES CALLED  SAXPY */
/* ***END PROLOGUE  SHEQR */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */

/*         Local Variables. */


/*         minimize(B-A*X,B-A*X).  First form Q*B. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --q;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	iq = (k - 1 << 1) + 1;
	c__ = q[iq];
	s = q[iq + 1];
	t1 = b[k];
	t2 = b[kp1];
	b[k] = c__ * t1 - s * t2;
	b[kp1] = s * t1 + c__ * t2;
/* L20: */
    }

/*         Now solve  R*X = Q*B. */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L40: */
    }
    return 0;
/* ------------- LAST LINE OF SHELS FOLLOWS ---------------------------- */
} /* shels_ */

/* DECK SHEQR */
/* Subroutine */ int sheqr_(real *a, integer *lda, integer *n, real *q, 
	integer *info, integer *ijob)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static real c__;
    static integer i__, j, k;
    static real s, t, t1, t2;
    static integer iq, km1, kp1, nm1;

/* ***BEGIN PROLOGUE  SHEQR */
/* ***DATE WRITTEN   871001   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SHEQR-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Generalized Minimum Residual */
/* ***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov */
/*           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/* ***PURPOSE  Internal routine for SGMRES. */
/* ***DESCRIPTION */
/*        This   routine  performs  a QR   decomposition  of an  upper */
/*        Hessenberg matrix A using Givens  rotations.  There  are two */
/*        options  available: 1)  Performing  a fresh decomposition 2) */
/*        updating the QR factors by adding a row and  a column to the */
/*        matrix A. */

/* *Usage: */
/*      INTEGER LDA, N, INFO, IJOB */
/*      REAL A(LDA,1), Q(1) */

/*      CALL SHEQR(A, LDA, N, Q, INFO, IJOB) */

/* *Arguments: */
/* A      :INOUT    Real A(LDA,N) */
/*         On input, the matrix to be decomposed. */
/*         On output, the upper triangular matrix R. */
/*         The factorization can be written Q*A = R, where */
/*         Q is a product of Givens rotations and R is upper */
/*         triangular. */
/* LDA    :IN       Integer */
/*         The leading dimension of the array A. */
/* N      :IN       Integer */
/*         A is an (N+1) by N Hessenberg matrix. */
/* IJOB   :IN       Integer */
/*         = 1     means that a fresh decomposition of the */
/*                 matrix A is desired. */
/*         .ge. 2  means that the current decomposition of A */
/*                 will be updated by the addition of a row */
/*                 and a column. */
/* Q      :OUT      Real Q(2*N) */
/*         The factors c and s of each Givens rotation used */
/*         in decomposing A. */
/* INFO   :OUT      Integer */
/*         = 0  normal value. */
/*         = K  if  A(K,K) .eq. 0.0 .  This is not an error */
/*           condition for this subroutine, but it does */
/*           indicate that SHELS will divide by zero */
/*           if called. */

/* *See Also: */
/*         SGMRES */

/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SHEQR */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */

/*         Local Variables. */


/* ***FIRST EXECUTABLE STATEMENT  SHEQR */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --q;

    /* Function Body */
    if (*ijob > 1) {
	goto L70;
    }
/*   ------------------------------------------------------------------- */
/*         A new facorization is desired. */
/*   ------------------------------------------------------------------- */
/*         QR decomposition without pivoting. */

    *info = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	km1 = k - 1;
	kp1 = k + 1;

/*           Compute K-th column of R. */
/*           First, multiply the K-th column of a by the previous */
/*           K-1 Givens rotations. */

	if (km1 < 1) {
	    goto L20;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	    i__ = (j - 1 << 1) + 1;
	    t1 = a[j + k * a_dim1];
	    t2 = a[j + 1 + k * a_dim1];
	    c__ = q[i__];
	    s = q[i__ + 1];
	    a[j + k * a_dim1] = c__ * t1 - s * t2;
	    a[j + 1 + k * a_dim1] = s * t1 + c__ * t2;
/* L10: */
	}

/*         Compute Givens components C and S. */

L20:
	iq = (km1 << 1) + 1;
	t1 = a[k + k * a_dim1];
	t2 = a[kp1 + k * a_dim1];
	if (t2 == 0.f) {
	    c__ = 1.f;
	    s = 0.f;
	} else if (dabs(t2) >= dabs(t1)) {
	    t = t1 / t2;
	    s = -1.f / sqrt(t * t + 1.f);
	    c__ = -s * t;
	} else {
	    t = t2 / t1;
	    c__ = 1.f / sqrt(t * t + 1.f);
	    s = -c__ * t;
	}
	q[iq] = c__;
	q[iq + 1] = s;
	a[k + k * a_dim1] = c__ * t1 - s * t2;
	if (a[k + k * a_dim1] == 0.f) {
	    *info = k;
	}
/* L60: */
    }
    return 0;
/*   ------------------------------------------------------------------- */
/*         The old factorization of a will be updated.  A row and a */
/*         column has been added to the matrix A.  N by N-1 is now */
/*         the old size of the matrix. */
/*   ------------------------------------------------------------------- */
L70:
    nm1 = *n - 1;
/*   ------------------------------------------------------------------- */
/*         Multiply the new column by the N previous Givens rotations. */
/*   ------------------------------------------------------------------- */
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	i__ = (k - 1 << 1) + 1;
	t1 = a[k + *n * a_dim1];
	t2 = a[k + 1 + *n * a_dim1];
	c__ = q[i__];
	s = q[i__ + 1];
	a[k + *n * a_dim1] = c__ * t1 - s * t2;
	a[k + 1 + *n * a_dim1] = s * t1 + c__ * t2;
/* L100: */
    }
/*   ------------------------------------------------------------------- */
/*         Complete update of decomposition by forming last Givens */
/*         rotation, and multiplying it times the column */
/*         vector(A(N,N),A(NP1,N)). */
/*   ------------------------------------------------------------------- */
    *info = 0;
    t1 = a[*n + *n * a_dim1];
    t2 = a[*n + 1 + *n * a_dim1];
    if (t2 == 0.f) {
	c__ = 1.f;
	s = 0.f;
    } else if (dabs(t2) >= dabs(t1)) {
	t = t1 / t2;
	s = -1.f / sqrt(t * t + 1.f);
	c__ = -s * t;
    } else {
	t = t2 / t1;
	c__ = 1.f / sqrt(t * t + 1.f);
	s = -c__ * t;
    }
    iq = (*n << 1) - 1;
    q[iq] = c__;
    q[iq + 1] = s;
    a[*n + *n * a_dim1] = c__ * t1 - s * t2;
    if (a[*n + *n * a_dim1] == 0.f) {
	*info = *n;
    }
    return 0;
/* ------------- LAST LINE OF SHEQR FOLLOWS ---------------------------- */
} /* sheqr_ */

/* DECK SORTH */
/* Subroutine */ int sorth_(real *vnew, real *v, real *hes, integer *n, 
	integer *ll, integer *ldhes, integer *kmp, real *snormw)
{
    /* System generated locals */
    integer v_dim1, v_offset, hes_dim1, hes_offset, i__1, i__2;
    real r__1, r__2, r__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, i0;
    static real arg, tem;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    static real vnrm;
    extern doublereal snrm2_(integer *, real *, integer *);
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);
    static real sumdsq;

/* ***BEGIN PROLOGUE  SORTH */
/* ***DATE WRITTEN   871001   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SORTH-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Generalized Minimum Residual */
/* ***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov */
/*           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/* ***PURPOSE  Internal routine for SGMRES. */
/* ***DESCRIPTION */
/*        This routine  orthogonalizes  the  vector  VNEW  against the */
/*        previous KMP  vectors in the   V array.  It uses  a modified */
/*        gram-schmidt   orthogonalization procedure with  conditional */
/*        reorthogonalization. */

/* *Usage: */
/*      INTEGER N, LL, LDHES, KMP */
/*      REAL VNEW, V, HES, SNORMW */
/*      DIMENSION VNEW(1), V(N,1), HES(LDHES,1) */

/*      CALL SORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW) */

/* *Arguments: */
/* VNEW   :INOUT    Real VNEW(N) */
/*         On input, the vector of length n containing a scaled */
/*         product of the jacobian and the vector v(*,ll). */
/*         On output, the new vector orthogonal to v(*,i0) to v(*,ll), */
/*         where i0 = max(1, ll-kmp+1). */
/* V      :IN       Real V(N,1) */
/*         The n x ll array containing the previous ll */
/*         orthogonal vectors v(*,1) to v(*,ll). */
/* HES    :INOUT    Real HES(LDHES,1) */
/*         On input, an LL x LL upper hessenberg matrix containing, */
/*         in HES(I,K), K.lt.LL, the scaled inner products of */
/*         A*V(*,K) and V(*,i). */
/*         On return, column LL of HES is filled in with */
/*         the scaled inner products of A*V(*,LL) and V(*,i). */
/* LDHES  :IN       Integer */
/*         The leading dimension of the HES array. */
/* N      :IN       Integer */
/*         The order of the matrix A, and the length of VNEW. */
/* LL     :IN       Integer */
/*         The current order of the matrix HES. */
/* KMP    :IN       Integer */
/*         The number of previous vectors the new vector VNEW */
/*         must be made orthogonal to (KMP .le. MAXL). */
/* SNORMW :OUT      REAL */
/*         Scalar containing the l-2 norm of VNEW. */

/* *See Also: */
/*         SGMRES */

/* ***ROUTINES CALLED  SAXPY */
/* ***END PROLOGUE  SORTH */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */

/*         Internal variables. */


/*         Get norm of unaltered VNEW for later use. */
/* ***FIRST EXECUTABLE STATEMENT  SORTH */
    /* Parameter adjustments */
    --vnew;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    hes_dim1 = *ldhes;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;

    /* Function Body */
    vnrm = snrm2_(n, &vnew[1], &c__1);
/*   ------------------------------------------------------------------- */
/*         Perform the modified gram-schmidt procedure on VNEW =A*V(LL). */
/*         Scaled inner products give new column of HES. */
/*         Projections of earlier vectors are subtracted from VNEW. */
/*   ------------------------------------------------------------------- */
/* Computing MAX */
    i__1 = 1, i__2 = *ll - *kmp + 1;
    i0 = max(i__1,i__2);
    i__1 = *ll;
    for (i__ = i0; i__ <= i__1; ++i__) {
	hes[i__ + *ll * hes_dim1] = sdot_(n, &v[i__ * v_dim1 + 1], &c__1, &
		vnew[1], &c__1);
	tem = -hes[i__ + *ll * hes_dim1];
	saxpy_(n, &tem, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
/* L10: */
    }
/*   ------------------------------------------------------------------- */
/*         Compute SNORMW = norm of VNEW.  If VNEW is small compared */
/*         to its input value (in norm), then reorthogonalize VNEW to */
/*         V(*,1) through V(*,LL).  Correct if relative correction */
/*         exceeds 1000*(unit roundoff).  Finally, correct SNORMW using */
/*         the dot products involved. */
/*   ------------------------------------------------------------------- */
    *snormw = snrm2_(n, &vnew[1], &c__1);
    if (vnrm + *snormw * .001f != vnrm) {
	return 0;
    }
    sumdsq = 0.f;
    i__1 = *ll;
    for (i__ = i0; i__ <= i__1; ++i__) {
	tem = -sdot_(n, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
	if (hes[i__ + *ll * hes_dim1] + tem * .001f == hes[i__ + *ll * 
		hes_dim1]) {
	    goto L30;
	}
	hes[i__ + *ll * hes_dim1] -= tem;
	saxpy_(n, &tem, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
/* Computing 2nd power */
	r__1 = tem;
	sumdsq += r__1 * r__1;
L30:
	;
    }
    if (sumdsq == 0.f) {
	return 0;
    }
/* Computing MAX */
/* Computing 2nd power */
    r__3 = *snormw;
    r__1 = 0.f, r__2 = r__3 * r__3 - sumdsq;
    arg = dmax(r__1,r__2);
    *snormw = sqrt(arg);

    return 0;
/* ------------- LAST LINE OF SORTH FOLLOWS ---------------------------- */
} /* sorth_ */

/* DECK SPIGMR */
/* Subroutine */ int spigmr_(integer *n, real *r0, real *sr, real *sz, 
	integer *jscal, integer *maxl, integer *maxlp1, integer *kmp, integer 
	*nrsts, integer *jpre, S_fp matvec, S_fp msolve, integer *nmsl, real *
	z__, real *v, real *hes, real *q, integer *lgmr, real *rpar, integer *
	ipar, real *wk, real *dl, real *rhol, integer *nrmax, real *b, real *
	bnrm, real *x, real *xl, integer *itol, real *tol, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, integer *iunit, 
	integer *iflag, real *err)
{
    /* System generated locals */
    integer v_dim1, v_offset, hes_dim1, hes_offset, i__1, i__2, i__3;
    real r__1;

    /* Local variables */
    static real c__;
    static integer i__, j, k;
    static real s;
    static integer i2, ll, ip1;
    static real tem, rho;
    static integer llp1, info;
    static real prod;
    static integer iter;
    static real r0nrm;
    extern doublereal snrm2_(integer *, real *, integer *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static real dlnrm;
    extern /* Subroutine */ int shels_(real *, integer *, integer *, real *, 
	    real *);
    static integer itmax;
    extern /* Subroutine */ int sheqr_(real *, integer *, integer *, real *, 
	    integer *, integer *), scopy_(integer *, real *, integer *, real *
	    , integer *), sorth_(real *, real *, real *, integer *, integer *,
	     integer *, integer *, real *), saxpy_(integer *, real *, real *, 
	    integer *, real *, integer *), srlcal_(integer *, integer *, 
	    integer *, integer *, real *, real *, real *, real *, real *, 
	    real *);
    extern integer issgmr_(integer *, real *, real *, real *, integer *, 
	    integer *, integer *, real *, integer *, S_fp, integer *, integer 
	    *, real *, integer *, integer *, real *, integer *, real *, real *
	    , real *, real *, integer *, real *, real *, real *, real *, 
	    integer *, integer *, integer *, integer *, integer *, real *, 
	    real *, real *, real *, real *, real *, integer *);
    static real snormw;

/* ***BEGIN PROLOGUE  SPIGMR */
/* ***DATE WRITTEN   871001   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SPIGMR-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Generalized Minimum Residual */
/* ***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov */
/*           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/* ***PURPOSE  Internal routine for SGMRES. */
/* ***DESCRIPTION */
/*         This routine solves the linear system A * Z = R0 using a */
/*         scaled preconditioned version of the generalized minimum */
/*         residual method.  An initial guess of Z = 0 is assumed. */

/* *Usage: */
/*      EXTERNAL MATVEC, MSOLVE */
/*      INTEGER N,MAXL,MAXLP1,KMP,JPRE,NMSL,LGMR,IPAR,IFLAG,JSCAL,NRSTS */
/*      INTEGER NRMAX,ITOL,NELT,ISYM */
/*      REAL R0,SR,SZ,Z,V,HES,Q,RPAR,WK,DL,RHOL,BNRM,TOL,A,B,X */
/*      DIMENSION R0(1), SR(1), SZ(1), Z(1), V(N,1), */
/*     $     HES(MAXLP1,1), Q(1), RPAR(1), IPAR(1), WK(1), DL(1), */
/*     $     IA(NELT), JA(NELT), A(NELT), B(1), X(1), XL(1) */

/*      CALL SPIGMR(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, */
/*     $     NRSTS, JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR, */
/*     $     RPAR, IPAR, WK, DL, RHOL, NRMAX, B, BNRM, X, XL, */
/*     $     ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR) */

/* *Arguments: */
/* R0     :IN       Real R0(N) */
/*         R0 = the right hand side of the system A*Z = R0. */
/*         R0 is also used as work space when computing */
/*         the final approximation. */
/*         (R0 is the same as V(*,MAXL+1) in the call to SPIGMR.) */
/* SR     :IN       Real SR(N) */
/*         SR is a vector of length N containing the nonzero */
/*         elements of the diagonal scaling matrix for R0. */
/* SZ     :IN       Real SZ(N) */
/*         SZ is a vector of length N containing the nonzero */
/*         elements of the diagonal scaling matrix for Z. */
/* JSCAL  :IN       Integer */
/*         A flag indicating whether arrays SR and SZ are used. */
/*         JSCAL=0 means SR and SZ are not used and the */
/*                 algorithm will perform as if all */
/*                 SR(i) = 1 and SZ(i) = 1. */
/*         JSCAL=1 means only SZ is used, and the algorithm */
/*                 performs as if all SR(i) = 1. */
/*         JSCAL=2 means only SR is used, and the algorithm */
/*                 performs as if all SZ(i) = 1. */
/*         JSCAL=3 means both SR and SZ are used. */
/* N      :IN       Integer */
/*         The order of the matrix A, and the lengths */
/*         of the vectors SR, SZ, R0 and Z. */
/* MAXL   :IN       Integer */
/*         The maximum allowable order of the matrix H. */
/* MAXLP1 :IN       Integer */
/*         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES. */
/* KMP    :IN       Integer */
/*         The number of previous vectors the new vector VNEW */
/*         must be made orthogonal to.  (KMP .le. MAXL) */
/* NRSTS  :IN       Integer */
/*         Counter for the number of restarts on the current */
/*         call to SGMRES.  If NRSTS .gt. 0, then the residual */
/*         R0 is already scaled, and so scaling of it is */
/*         not necessary. */
/* JPRE   :IN       Integer */
/*         Preconditioner type flag. */
/* WK     :IN       Real WK(N) */
/*         A real work array of length N used by routine MATVEC */
/*         and MSOLVE. */
/* DL     :INOUT    Real DL(N) */
/*         On input, a real work array of length N used for calculation */
/*         of the residual norm RHO when the method is incomplete */
/*         (KMP.lt.MAXL), and/or when using restarting. */
/*         On output, the scaled residual vector RL.  It is only loaded */
/*         when performing restarts of the Krylov iteration. */
/* NRMAX  :IN       Integer */
/*         The maximum number of restarts of the Krylov iteration. */
/*         NRMAX .gt. 0 means restarting is active, while */
/*         NRMAX = 0 means restarting is not being used. */
/* B      :IN       Real B(N) */
/*         The right hand side of the linear system A*X = B. */
/* BNRM   :IN       Real */
/*         The scaled norm of b. */
/* X      :IN       Real X(N) */
/*         The current approximate solution as of the last */
/*         restart. */
/* XL     :IN       Real XL(N) */
/*         An array of length N used to hold the approximate */
/*         solution X(L) when ITOL=11. */
/* ITOL   :IN       Integer */
/*         A flag to indicate the type of convergence criterion */
/*         used.  see the driver for its description. */
/* TOL    :IN       Real */
/*         The tolerance on residuals R0-A*Z in scaled norm. */
/* NELT   :IN       Integer */
/*         The length of arrays IA, JA and A. */
/* IA     :IN       Integer IA(NELT) */
/*         An integer array of length NELT containing matrix data. */
/*         It is passed directly to the MATVEC and MSOLVE routines. */
/* JA     :IN       Integer JA(NELT) */
/*         An integer array of length NELT containing matrix data. */
/*         It is passed directly to the MATVEC and MSOLVE routines. */
/* A      :IN       Real A(NELT) */
/*         A real array of length NELT containing matrix data. */
/*         It is passed directly to the MATVEC and MSOLVE routines. */
/* ISYM   :IN       Integer */
/*         A flag to indicate symmetric matrix storage. */
/*         If ISYM=0, all nonzero entries of the matrix are */
/*         stored.  If ISYM=1, the matrix is symmetric and */
/*         only the upper or lower triangular part is stored. */
/* IUNIT  :IN       Integer */
/*         The i/o unit number for writing intermediate residual */
/*         norm values. */
/* Z      :OUT      Real Z(N) */
/*         The final computed approximation to the solution */
/*         of the system A*Z = R0. */
/* LGMR   :OUT      Integer */
/*         The number of iterations performed and */
/*         the current order of the upper hessenberg */
/*         matrix HES. */
/* RPAR   :IN       Real RPAR(*) */
/*         Real work space passed directly to the MSOLVE routine. */
/* IPAR   :IN       Integer IPAR(*) */
/*         Integer work space passed directly to the MSOLVE */
/*         routine. */
/* NMSL   :OUT      Integer */
/*         The number of calls to MSOLVE. */
/* V      :OUT      Real V(N,MAXLP1) */
/*         The N by (LGMR+1) array containing the LGMR */
/*         orthogonal vectors V(*,1) to V(*,LGMR). */
/* HES    :OUT      Real HES(MAXLP1,MAXL) */
/*         The upper triangular factor of the QR decomposition */
/*         of the (LGMR+1) by LGMR upper Hessenberg matrix whose */
/*         entries are the scaled inner-products of A*V(*,I) */
/*         and V(*,K). */
/* Q      :OUT      Real Q(2*MAXL) */
/*         A real array of length 2*MAXL containing the components */
/*         of the Givens rotations used in the QR decomposition */
/*         of HES.  It is loaded in SHEQR and used in SHELS. */
/* RHOL   :OUT      Real */
/*         A real scalar containing the norm of the final residual. */
/* IFLAG  :OUT      Integer */
/*         An integer error flag.. */
/*         0 means convergence in LGMR iterations, LGMR.le.MAXL. */
/*         1 means the convergence test did not pass in MAXL */
/*           iterations, but the residual norm is .lt. norm(R0), */
/*           and so Z is computed. */
/*         2 means the convergence test did not pass in MAXL */
/*           iterations, residual .ge. norm(R0), and Z = 0. */
/* ERR    :OUT      Real. */
/*         Error estimate of error in final approximate solution, as */
/*         defined by ITOL. */

/* *See Also: */
/*         SGMRES */

/* ***ROUTINES CALLED  ISSGMR, MATVEC, MSOLVE, SORTH, SRLCAL, SHELS, */
/*                    SHEQR, SXLCAL, SAXPY, SCOPY, SSCAL, */
/* ***END PROLOGUE  SPIGMR */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */

/*         Local variables. */


/*         Zero out the z array. */
/* ***FIRST EXECUTABLE STATEMENT  SPIGMR */
    /* Parameter adjustments */
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --r0;
    --sr;
    --sz;
    hes_dim1 = *maxlp1;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;
    --z__;
    --q;
    --rpar;
    --ipar;
    --wk;
    --dl;
    --b;
    --x;
    --xl;
    --a;
    --ja;
    --ia;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__] = 0.f;
/* L5: */
    }

    *iflag = 0;
    *lgmr = 0;
    *nmsl = 0;
/*         Load ITMAX, the maximum number of iterations. */
    itmax = (*nrmax + 1) * *maxl;
/*   ------------------------------------------------------------------- */
/*         The initial residual is the vector R0. */
/*         Apply left precon. if JPRE < 0 and this is not a restart. */
/*         Apply scaling to R0 if JSCAL = 2 or 3. */
/*   ------------------------------------------------------------------- */
    if (*jpre < 0 && *nrsts == 0) {
	scopy_(n, &r0[1], &c__1, &wk[1], &c__1);
	(*msolve)(n, &wk[1], &r0[1], nelt, &ia[1], &ja[1], &a[1], isym, &rpar[
		1], &ipar[1]);
	++(*nmsl);
    }
    if ((*jscal == 2 || *jscal == 3) && *nrsts == 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    v[i__ + v_dim1] = r0[i__] * sr[i__];
/* L10: */
	}
    } else {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    v[i__ + v_dim1] = r0[i__];
/* L20: */
	}
    }
    r0nrm = snrm2_(n, &v[v_offset], &c__1);
    iter = *nrsts * *maxl;

/*         Call stopping routine ISSGMR. */

    if (issgmr_(n, &b[1], &x[1], &xl[1], nelt, &ia[1], &ja[1], &a[1], isym, (
	    S_fp)msolve, nmsl, itol, tol, &itmax, &iter, err, iunit, &v[
	    v_dim1 + 1], &z__[1], &wk[1], &rpar[1], &ipar[1], &r0nrm, bnrm, &
	    sr[1], &sz[1], jscal, kmp, lgmr, maxl, maxlp1, &v[v_offset], &q[1]
	    , &snormw, &prod, &r0nrm, &hes[hes_offset], jpre) != 0) {
	return 0;
    }
    tem = 1.f / r0nrm;
    sscal_(n, &tem, &v[v_dim1 + 1], &c__1);

/*         Zero out the HES array. */

    i__1 = *maxl;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *maxlp1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    hes[i__ + j * hes_dim1] = 0.f;
/* L40: */
	}
/* L50: */
    }
/*   ------------------------------------------------------------------- */
/*         main loop to compute the vectors V(*,2) to V(*,MAXL). */
/*         The running product PROD is needed for the convergence test. */
/*   ------------------------------------------------------------------- */
    prod = 1.f;
    i__1 = *maxl;
    for (ll = 1; ll <= i__1; ++ll) {
	*lgmr = ll;
/*   ------------------------------------------------------------------- */
/*        Unscale  the  current V(LL)  and store  in WK.  Call routine */
/*        msolve    to   compute(M-inverse)*WK,   where    M   is  the */
/*        preconditioner matrix.  Save the answer in Z.   Call routine */
/*        MATVEC to compute  VNEW  = A*Z,  where  A is  the the system */
/*        matrix.  save the answer in  V(LL+1).  Scale V(LL+1).   Call */
/*        routine SORTH  to  orthogonalize the    new vector VNEW   = */
/*        V(*,LL+1).  Call routine SHEQR to update the factors of HES. */
/*   ------------------------------------------------------------------- */
	if (*jscal == 1 || *jscal == 3) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		wk[i__] = v[i__ + ll * v_dim1] / sz[i__];
/* L60: */
	    }
	} else {
	    scopy_(n, &v[ll * v_dim1 + 1], &c__1, &wk[1], &c__1);
	}
	if (*jpre > 0) {
	    (*msolve)(n, &wk[1], &z__[1], nelt, &ia[1], &ja[1], &a[1], isym, &
		    rpar[1], &ipar[1]);
	    ++(*nmsl);
	    (*matvec)(n, &z__[1], &v[(ll + 1) * v_dim1 + 1], nelt, &ia[1], &
		    ja[1], &a[1], isym);
	} else {
	    (*matvec)(n, &wk[1], &v[(ll + 1) * v_dim1 + 1], nelt, &ia[1], &ja[
		    1], &a[1], isym);
	}
	if (*jpre < 0) {
	    scopy_(n, &v[(ll + 1) * v_dim1 + 1], &c__1, &wk[1], &c__1);
	    (*msolve)(n, &wk[1], &v[(ll + 1) * v_dim1 + 1], nelt, &ia[1], &ja[
		    1], &a[1], isym, &rpar[1], &ipar[1]);
	    ++(*nmsl);
	}
	if (*jscal == 2 || *jscal == 3) {
	    i__2 = *n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		v[i__ + (ll + 1) * v_dim1] *= sr[i__];
/* L65: */
	    }
	}
	sorth_(&v[(ll + 1) * v_dim1 + 1], &v[v_offset], &hes[hes_offset], n, &
		ll, maxlp1, kmp, &snormw);
	hes[ll + 1 + ll * hes_dim1] = snormw;
	sheqr_(&hes[hes_offset], maxlp1, &ll, &q[1], &info, &ll);
	if (info == ll) {
	    goto L120;
	}
/*   ------------------------------------------------------------------- */
/*         Update RHO, the estimate of the norm of the residual R0-A*ZL. */
/*         If KMP <  MAXL, then the vectors V(*,1),...,V(*,LL+1) are not */
/*         necessarily orthogonal for LL > KMP.  The vector DL must then */
/*         be computed, and its norm used in the calculation of RHO. */
/*   ------------------------------------------------------------------- */
	prod *= q[ll * 2];
	rho = (r__1 = prod * r0nrm, dabs(r__1));
	if (ll > *kmp && *kmp < *maxl) {
	    if (ll == *kmp + 1) {
		scopy_(n, &v[v_dim1 + 1], &c__1, &dl[1], &c__1);
		i__2 = *kmp;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ip1 = i__ + 1;
		    i2 = i__ << 1;
		    s = q[i2];
		    c__ = q[i2 - 1];
		    i__3 = *n;
		    for (k = 1; k <= i__3; ++k) {
			dl[k] = s * dl[k] + c__ * v[k + ip1 * v_dim1];
/* L70: */
		    }
/* L75: */
		}
	    }
	    s = q[ll * 2];
	    c__ = q[(ll << 1) - 1] / snormw;
	    llp1 = ll + 1;
	    i__2 = *n;
	    for (k = 1; k <= i__2; ++k) {
		dl[k] = s * dl[k] + c__ * v[k + llp1 * v_dim1];
/* L80: */
	    }
	    dlnrm = snrm2_(n, &dl[1], &c__1);
	    rho *= dlnrm;
	}
	*rhol = rho;
/*   ------------------------------------------------------------------- */
/*         Test for convergence.  If passed, compute approximation ZL. */
/*         If failed and LL < MAXL, then continue iterating. */
/*   ------------------------------------------------------------------- */
	iter = *nrsts * *maxl + *lgmr;
	if (issgmr_(n, &b[1], &x[1], &xl[1], nelt, &ia[1], &ja[1], &a[1], 
		isym, (S_fp)msolve, nmsl, itol, tol, &itmax, &iter, err, 
		iunit, &dl[1], &z__[1], &wk[1], &rpar[1], &ipar[1], rhol, 
		bnrm, &sr[1], &sz[1], jscal, kmp, lgmr, maxl, maxlp1, &v[
		v_offset], &q[1], &snormw, &prod, &r0nrm, &hes[hes_offset], 
		jpre) != 0) {
	    goto L200;
	}
	if (ll == *maxl) {
	    goto L100;
	}
/*   ------------------------------------------------------------------- */
/*         Rescale so that the norm of V(1,LL+1) is one. */
/*   ------------------------------------------------------------------- */
	tem = 1.f / snormw;
	sscal_(n, &tem, &v[(ll + 1) * v_dim1 + 1], &c__1);
/* L90: */
    }
L100:
    if (rho < r0nrm) {
	goto L150;
    }
L120:
    *iflag = 2;

/*         Load approximate solution with zero. */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	z__[i__] = 0.f;
/* L130: */
    }
    return 0;
L150:
    *iflag = 1;

/*         Tolerance not met, but residual norm reduced. */

    if (*nrmax > 0) {

/*        If performing restarting (NRMAX > 0)  calculate the residual */
/*        vector RL and  store it in the DL  array.  If the incomplete */
/*        version is being used (KMP < MAXL) then DL has  already been */
/*        calculated up to a scaling factor.   Use SRLCAL to calculate */
/*        the scaled residual vector. */

	srlcal_(n, kmp, maxl, maxl, &v[v_offset], &q[1], &dl[1], &snormw, &
		prod, &r0nrm);
    }
/*   ------------------------------------------------------------------- */
/*         Compute the approximation ZL to the solution.  Since the */
/*         vector Z was used as work space, and the initial guess */
/*         of the linear iteration is zero, Z must be reset to zero. */
/*   ------------------------------------------------------------------- */
L200:
    ll = *lgmr;
    llp1 = ll + 1;
    i__1 = llp1;
    for (k = 1; k <= i__1; ++k) {
	r0[k] = 0.f;
/* L210: */
    }
    r0[1] = r0nrm;
    shels_(&hes[hes_offset], maxlp1, &ll, &q[1], &r0[1]);
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	z__[k] = 0.f;
/* L220: */
    }
    i__1 = ll;
    for (i__ = 1; i__ <= i__1; ++i__) {
	saxpy_(n, &r0[i__], &v[i__ * v_dim1 + 1], &c__1, &z__[1], &c__1);
/* L230: */
    }
    if (*jscal == 1 || *jscal == 3) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    z__[i__] /= sz[i__];
/* L240: */
	}
    }
    if (*jpre > 0) {
	scopy_(n, &z__[1], &c__1, &wk[1], &c__1);
	(*msolve)(n, &wk[1], &z__[1], nelt, &ia[1], &ja[1], &a[1], isym, &
		rpar[1], &ipar[1]);
	++(*nmsl);
    }
    return 0;
/* ------------- LAST LINE OF SPIGMR FOLLOWS ---------------------------- */
} /* spigmr_ */

/* DECK SRLCAL */
/* Subroutine */ int srlcal_(integer *n, integer *kmp, integer *ll, integer *
	maxl, real *v, real *q, real *rl, real *snormw, real *prod, real *
	r0nrm)
{
    /* System generated locals */
    integer v_dim1, v_offset, i__1, i__2;

    /* Local variables */
    static real c__;
    static integer i__, k;
    static real s;
    static integer i2, ip1;
    static real tem;
    static integer llm1, llp1;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *), 
	    scopy_(integer *, real *, integer *, real *, integer *);

/* ***BEGIN PROLOGUE  SRLCAL */
/* ***DATE WRITTEN   871001   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SRLCAL-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Generalized Minimum Residual */
/* ***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov */
/*           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/* ***PURPOSE  Internal routine for SGMRES. */
/* ***DESCRIPTION */
/*         This routine calculates the scaled residual RL from the */
/*         V(I)'s. */
/* *Usage: */
/*      INTEGER N, KMP, LL, MAXL */
/*      REAL V, Q, RL, SNORMW */
/*      DIMENSION V(N,1), Q(1), RL(N) */

/*      CALL SRLCAL(N, KMP, LL, MAXL, V, Q, RL, SNORMW, PROD, */
/*     $     R0NRM) */

/* *Arguments: */
/* N      :IN       Integer */
/*         The order of the matrix A, and the lengths */
/*         of the vectors SR, SZ, R0 and Z. */
/* KMP    :IN       Integer */
/*         The number of previous V vectors the new vector VNEW */
/*         must be made orthogonal to. (KMP .le. MAXL) */
/* LL     :IN       Integer */
/*         The current dimension of the Krylov subspace. */
/* MAXL   :IN       Integer */
/*         The maximum dimension of the Krylov subspace. */
/* Q      :IN       Real Q(2*MAXL) */
/*         A real array of length 2*MAXL containing the components */
/*         of the Givens rotations used in the QR decomposition */
/*         of HES.  It is loaded in SHEQR and used in SHELS. */
/* PROD   :IN       Real */
/*        The product s1*s2*...*sl = the product of the sines of the */
/*        givens rotations used in the QR factorization of */
/*        the hessenberg matrix HES. */
/* R0NRM  :IN       Real */
/*         The scaled norm of initial residual R0. */
/* RL     :OUT      Real RL(N) */
/*         The residual vector RL.  This is either SB*(B-A*XL) if */
/*         not preconditioning or preconditioning on the right, */
/*         or SB*(M-inverse)*(B-A*XL) if preconditioning on the */
/*         left. */

/* *See Also: */
/*         SGMRES */

/* ***ROUTINES CALLED  SCOPY, SSCAL */
/* ***END PROLOGUE  SRLCAL */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */

/*         Internal Variables. */


/* ***FIRST EXECUTABLE STATEMENT  SRLCAL */
    /* Parameter adjustments */
    --rl;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --q;

    /* Function Body */
    if (*kmp == *maxl) {

/*         calculate RL.  Start by copying V(*,1) into RL. */

	scopy_(n, &v[v_dim1 + 1], &c__1, &rl[1], &c__1);
	llm1 = *ll - 1;
	i__1 = llm1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ip1 = i__ + 1;
	    i2 = i__ << 1;
	    s = q[i2];
	    c__ = q[i2 - 1];
	    i__2 = *n;
	    for (k = 1; k <= i__2; ++k) {
		rl[k] = s * rl[k] + c__ * v[k + ip1 * v_dim1];
/* L10: */
	    }
/* L20: */
	}
	s = q[*ll * 2];
	c__ = q[(*ll << 1) - 1] / *snormw;
	llp1 = *ll + 1;
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    rl[k] = s * rl[k] + c__ * v[k + llp1 * v_dim1];
/* L30: */
	}
    }

/*         When KMP < MAXL, RL vector already partially calculated. */
/*         Scale RL by R0NRM*PROD to obtain the residual RL. */

    tem = *r0nrm * *prod;
    sscal_(n, &tem, &rl[1], &c__1);
    return 0;
/* ------------- LAST LINE OF SRLCAL FOLLOWS ---------------------------- */
} /* srlcal_ */

/* DECK SXLCAL */
/* Subroutine */ int sxlcal_(integer *n, integer *lgmr, real *x, real *xl, 
	real *zl, real *hes, integer *maxlp1, real *q, real *v, real *r0nrm, 
	real *wk, real *sz, integer *jscal, integer *jpre, S_fp msolve, 
	integer *nmsl, real *rpar, integer *ipar, integer *nelt, integer *ia, 
	integer *ja, real *a, integer *isym)
{
    /* System generated locals */
    integer hes_dim1, hes_offset, v_dim1, v_offset, i__1;

    /* Local variables */
    static integer i__, k, ll, llp1;
    extern /* Subroutine */ int shels_(real *, integer *, integer *, real *, 
	    real *), scopy_(integer *, real *, integer *, real *, integer *), 
	    saxpy_(integer *, real *, real *, integer *, real *, integer *);

/* ***BEGIN PROLOGUE  SXLCAL */
/* ***DATE WRITTEN   871001   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SXLCAL-S), */
/*             Non-Symmetric Linear system, Sparse, */
/*             Iterative Precondition, Generalized Minimum Residual */
/* ***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov */
/*           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/* ***PURPOSE  Internal routine for SGMRES. */
/* ***DESCRIPTION */
/*        This  routine computes the solution  XL,  the current SGMRES */
/*        iterate, given the  V(I)'s and  the  QR factorization of the */
/*        Hessenberg  matrix HES.   This routine  is  only called when */
/*        ITOL=11. */

/* *Usage: */
/*      EXTERNAL MSOLVE */
/*      INTEGER N, LGMR, MAXLP1, JSCAL, JPRE, IPAR, NMSL, NELT, IA, */
/*     $        JA, ISYM */
/*      REAL X, XL, ZL, HES, Q, V, R0NRM, WK, SZ, RPAR, A */
/*      DIMENSION X(N), XL(N), ZL(N), HES(MAXLP1,1), Q(1), V(N,1), WK(N), */
/*     $          SZ(1), RPAR(1), IPAR(1), IA(NELT), JA(NELT), A(NELT) */
/*      CALL SXLCAL(N, LGMR, X, XL, ZL, HES, MAXLP1, Q, V, R0NRM, */
/*     $     WK, SZ, JSCAL, JPRE, MSOLVE, NMSL, RPAR, IPAR, */
/*     $     NELT, IA, JA, A, ISYM) */

/* *Arguments: */
/* N      :IN       Integer */
/*         The order of the matrix A, and the lengths */
/*         of the vectors SR, SZ, R0 and Z. */
/* LGMR   :IN       Integer */
/*         The number of iterations performed and */
/*         the current order of the upper Hessenberg */
/*         matrix HES. */
/* X      :IN       Real X(N) */
/*         The current approximate solution as of the last restart. */
/* ZL     :IN       Real ZL(N) */
/*         An array of length N used to hold the approximate */
/*         solution Z(L). */
/* SZ     :IN       Real SZ(N) */
/*         A vector of length N containing the nonzero */
/*         elements of the diagonal scaling matrix for Z. */
/* JSCAL  :IN       Integer */
/*         A flag indicating whether arrays SR and SZ are used. */
/*         JSCAL=0 means SR and SZ are not used and the */
/*                 algorithm will perform as if all */
/*                 SR(i) = 1 and SZ(i) = 1. */
/*         JSCAL=1 means only SZ is used, and the algorithm */
/*                 performs as if all SR(i) = 1. */
/*         JSCAL=2 means only SR is used, and the algorithm */
/*                 performs as if all SZ(i) = 1. */
/*         JSCAL=3 means both SR and SZ are used. */
/* MAXLP1 :IN       Integer */
/*         MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES. */
/*         MAXL is the maximum allowable order of the matrix HES. */
/* JPRE   :IN       Integer */
/*         The preconditioner type flag. */
/* WK     :IN       Real WK(N) */
/*         A real work array of length N. */
/* NMSL   :IN       Integer */
/*         The number of calls to MSOLVE. */
/* V      :IN       Real V(N,MAXLP1) */
/*         The N by(LGMR+1) array containing the LGMR */
/*         orthogonal vectors V(*,1) to V(*,LGMR). */
/* HES    :IN       Real HES(MAXLP1,MAXL) */
/*         The upper triangular factor of the QR decomposition */
/*         of the (LGMR+1) by LGMR upper Hessenberg matrix whose */
/*         entries are the scaled inner-products of A*V(*,i) and V(*,k). */
/* Q      :IN       Real Q(2*MAXL) */
/*         A real array of length 2*MAXL containing the components */
/*         of the givens rotations used in the QR decomposition */
/*         of HES.  It is loaded in SHEQR. */
/* R0NRM  :IN       Real */
/*         The scaled norm of the initial residual for the */
/*         current call to SPIGMR. */
/* RPAR   :IN       Real RPAR(*) */
/*         Real work space passed directly to the MSOLVE routine. */
/* IPAR   :IN       Integer IPAR(*) */
/*         Integer work space passed directly to the MSOLVE */
/*         routine. */
/* NELT   :IN       Integer */
/*         The length of arrays IA, JA and A. */
/* IA     :IN       Integer IA(NELT) */
/*         An integer array of length NELT containing matrix data. */
/*         It is passed directly to the MATVEC and MSOLVE routines. */
/* JA     :IN       Integer JA(NELT) */
/*         An integer array of length NELT containing matrix data. */
/*         It is passed directly to the MATVEC and MSOLVE routines. */
/* A      :IN       Real A(NELT) */
/*         A real array of length NELT containing matrix data. */
/*         It is passed directly to the MATVEC and MSOLVE routines. */
/* ISYM   :IN       Integer */
/*         A flag to indicate symmetric matrix storage. */
/*         If ISYM=0, all nonzero entries of the matrix are */
/*         stored.  If ISYM=1, the matrix is symmetric and */
/*         only the upper or lower triangular part is stored. */
/* XL     :OUT      Real XL(N) */
/*         An array of length N used to hold the approximate */
/*         solution X(L). */
/*         Warning: XL and ZL are the same array in the calling routine. */

/* *See Also: */
/*         SGMRES */

/* ***ROUTINES CALLED  MSOLVE, SHELS, SAXPY, SCOPY, SSCAL */
/* ***END PROLOGUE  SXLCAL */
/*         The following is for optimized compilation on LLNL/LTSS Crays. */
/* LLL. OPTIMIZE */

/*         Internal variables. */


/* ***FIRST EXECUTABLE STATEMENT  SXLCAL */
    /* Parameter adjustments */
    --wk;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --zl;
    --xl;
    --x;
    hes_dim1 = *maxlp1;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;
    --q;
    --sz;
    --rpar;
    --ipar;
    --a;
    --ja;
    --ia;

    /* Function Body */
    ll = *lgmr;
    llp1 = ll + 1;
    i__1 = llp1;
    for (k = 1; k <= i__1; ++k) {
	wk[k] = 0.f;
/* L10: */
    }
    wk[1] = *r0nrm;
    shels_(&hes[hes_offset], maxlp1, &ll, &q[1], &wk[1]);
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	zl[k] = 0.f;
/* L20: */
    }
    i__1 = ll;
    for (i__ = 1; i__ <= i__1; ++i__) {
	saxpy_(n, &wk[i__], &v[i__ * v_dim1 + 1], &c__1, &zl[1], &c__1);
/* L30: */
    }
    if (*jscal == 1 || *jscal == 3) {
	i__1 = *n;
	for (k = 1; k <= i__1; ++k) {
	    zl[k] /= sz[k];
/* L40: */
	}
    }
    if (*jpre > 0) {
	scopy_(n, &zl[1], &c__1, &wk[1], &c__1);
	(*msolve)(n, &wk[1], &zl[1], nelt, &ia[1], &ja[1], &a[1], isym, &rpar[
		1], &ipar[1]);
	++(*nmsl);
    }
/*         calculate XL from X and ZL. */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	xl[k] = x[k] + zl[k];
/* L50: */
    }
    return 0;
/* ------------- LAST LINE OF SXLCAL FOLLOWS ---------------------------- */
} /* sxlcal_ */

/* DECK ISSGMR */
integer issgmr_(integer *n, real *b, real *x, real *xl, integer *nelt, 
	integer *ia, integer *ja, real *a, integer *isym, S_fp msolve, 
	integer *nmsl, integer *itol, real *tol, integer *itmax, integer *
	iter, real *err, integer *iunit, real *r__, real *z__, real *dz, real 
	*rwork, integer *iwork, real *rnrm, real *bnrm, real *sb, real *sx, 
	integer *jscal, integer *kmp, integer *lgmr, integer *maxl, integer *
	maxlp1, real *v, real *q, real *snormw, real *prod, real *r0nrm, real 
	*hes, integer *jpre)
{
    /* Format strings */
    static char fmt_1020[] = "(1x,\002 ITER = \002,i5,\002 IELMAX = \002,i5"
	    ",\002 |R(IELMAX)/X(IELMAX)| = \002,e12.5)";
    static char fmt_1000[] = "(\002 Generalized Minimum Residual(\002,i3,i3"
	    ",\002) for \002,\002N, ITOL = \002,i5,i5,/\002 ITER\002,\002   N"
	    "atral Err Est\002,\002   Error Estimate\002)";
    static char fmt_1010[] = "(1x,i4,1x,e16.7,1x,e16.7)";

    /* System generated locals */
    integer v_dim1, v_offset, hes_dim1, hes_offset, ret_val, i__1;
    real r__1, r__2, r__3;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static real tem, rat, fuzz;
    extern doublereal snrm2_(integer *, real *, integer *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    static real dxnrm;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *);
    extern doublereal r1mach_(integer *);
    static integer ielmax;
    extern /* Subroutine */ int srlcal_(integer *, integer *, integer *, 
	    integer *, real *, real *, real *, real *, real *, real *), 
	    sxlcal_(integer *, integer *, real *, real *, real *, real *, 
	    integer *, real *, real *, real *, real *, real *, integer *, 
	    integer *, S_fp, integer *, real *, integer *, integer *, integer 
	    *, integer *, real *, integer *);
    static real ratmax, solnrm;

    /* Fortran I/O blocks */
    static cilist io___115 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___118 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___119 = { 0, 0, 0, fmt_1010, 0 };


/* ***BEGIN PROLOGUE ISSGMR */
/* ***DATE WRITTEN   871211  (YYMMDD) */
/* ***REVISION DATE  881213  (YYMMDD) */
/* ***CATEGORY NO.  D2A4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=INTEGER(ISSGMR-I) */
/*             Linear system, Sparse, Stop Test, GMRES */
/* ***AUTHOR  Brown, Peter,    (LLNL), brown@lll-crg.llnl.gov */
/*           Hindmarsh, Alan, (LLNL), alanh@lll-crg.llnl.gov */
/*           Seager, Mark K., (LLNL), seager@lll-crg.llnl.gov */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/* ***PURPOSE Generalized Minimum Residual Stop Test. */
/*         This routine calculates the stop test for the Generalized */
/*         Minimum RESidual (GMRES) iteration scheme.   It returns a */
/*         nonzero if  the  error  estimate (the  type  of  which is */
/*         determined  by   ITOL)  is  less  than the user specified */
/*         tolerence TOL. */
/* ***DESCRIPTION */
/* *Usage: */
/*      INTEGER KMP, LGMR, MAXL, MAXLP1, JPRE, NMSL */
/*      REAL DXNRM, RNRM, R0NRM, SNORMW, SOLNRM, PROD */
/*      DIMENSION B(1), X(1), IA(1), JA(1), A(1), R(1), Z(1), DZ(1) */
/*      DIMENSION RWORK(1), IWORK(1), SB(1), SX(1), Q(1), V(N,1) */
/*      DIMENSION HES(MAXLP1,MAXL), XL(1) */
/*      EXTERNAL MSOLVE */

/*      IF (ISSGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE, */
/*     $     NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ, */
/*     $     RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL, */
/*     $     KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM, */
/*     $     HES, JPRE) .NE. 0) THEN ITERATION DONE */

/* *Arguments: */
/* N      :IN       Integer. */
/*         Order of the Matrix. */
/* B      :IN       Real B(N). */
/*         Right-hand-side vector. */
/* X      :IN       Real X(N). */
/*         Approximate solution vector as of the last restart. */
/* XL     :OUT      Real XL(N) */
/*         An array of length N used to hold the approximate */
/*         solution as of the current iteration.  Only computed by */
/*         this routine when ITOL=11. */
/* NELT   :IN       Integer. */
/*         Number of Non-Zeros stored in A. */
/* IA     :IN       Integer IA(NELT). */
/* JA     :IN       Integer JA(NELT). */
/* A      :IN       Real A(NELT). */
/*         These arrays contain the matrix data structure for A. */
/*         It could take any form.  See "Description", in the SGMRES, */
/*         SSLUGM and SSDGMR routines for more late breaking details... */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the upper */
/*         or lower triangle of the matrix is stored. */
/* MSOLVE :EXT      External. */
/*         Name of a routine which solves a linear system Mz = r for  z */
/*         given r with the preconditioning matrix M (M is supplied via */
/*         RWORK and IWORK arrays.  The name of the MSOLVE routine must */
/*         be declared external in the calling program.  The calling */
/*         sequence to MSLOVE is: */
/*             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK) */
/*         Where N is the number of unknowns, R is the right-hand side */
/*         vector, and z is the solution upon return.  RWORK is a real */
/*         array that can be used to pass necessary preconditioning */
/*         information and/or workspace to MSOLVE.  IWORK is an integer */
/*         work array for the same purpose as RWORK. */
/* NMSL   :INOUT    Integer. */
/*         A counter for the number of calls to MSOLVE. */
/* ITOL   :IN       Integer. */
/*         Flag to indicate the type of convergence criterion used. */
/*         ITOL=0  Means the  iteration stops when the test described */
/*                 below on  the  residual RL  is satisfied.  This is */
/*                 the  "Natural Stopping Criteria" for this routine. */
/*                 Other values  of   ITOL  cause  extra,   otherwise */
/*                 unnecessary, computation per iteration and     are */
/*                 therefore  much less  efficient.  See  ISSGMR (the */
/*                 stop test routine) for more information. */
/*         ITOL=1  Means   the  iteration stops   when the first test */
/*                 described below on  the residual RL  is satisfied, */
/*                 and there  is either right  or  no preconditioning */
/*                 being used. */
/*         ITOL=2  Implies     that   the  user    is   using    left */
/*                 preconditioning, and the second stopping criterion */
/*                 below is used. */
/*         ITOL=3  Means the  iteration stops   when  the  third test */
/*                 described below on Minv*Residual is satisfied, and */
/*                 there is either left  or no  preconditioning begin */
/*                 used. */
/*         ITOL=11 is    often  useful  for   checking  and comparing */
/*                 different routines.  For this case, the  user must */
/*                 supply  the  "exact" solution or  a  very accurate */
/*                 approximation (one with  an  error much less  than */
/*                 TOL) through a common block, */
/*                     COMMON /SOLBLK/ SOLN(1) */
/*                 if ITOL=11, iteration stops when the 2-norm of the */
/*                 difference between the iterative approximation and */
/*                 the user-supplied solution  divided by the  2-norm */
/*                 of the  user-supplied solution  is  less than TOL. */
/*                 Note that this requires  the  user to  set up  the */
/*                 "COMMON     /SOLBLK/ SOLN(LENGTH)"  in the calling */
/*                 routine.  The routine with this declaration should */
/*                 be loaded before the stop test so that the correct */
/*                 length is used by  the loader.  This procedure  is */
/*                 not standard Fortran and may not work correctly on */
/*                 your   system (although  it  has  worked  on every */
/*                 system the authors have tried).  If ITOL is not 11 */
/*                 then this common block is indeed standard Fortran. */
/* TOL    :IN       Real. */
/*         Convergence criterion, as described above. */
/* ITMAX  :IN       Integer. */
/*         Maximum number of iterations. */
/* ITER   :IN       Integer. */
/*         The iteration for which to check for convergence. */
/* ERR    :OUT      Real. */
/*         Error estimate of error in final approximate solution, as */
/*         defined by ITOL.  Letting norm() denote the Euclidean */
/*         norm, ERR is defined as follows.. */

/*         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B), */
/*                               for right or no preconditioning, and */
/*                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/ */
/*                                norm(SB*(M-inverse)*B), */
/*                               for left preconditioning. */
/*         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B), */
/*                               since right or no preconditioning */
/*                               being used. */
/*         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/ */
/*                                norm(SB*(M-inverse)*B), */
/*                               since left preconditioning is being */
/*                               used. */
/*         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)| */
/*                               i=1,n */
/*         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN). */
/* IUNIT  :IN       Integer. */
/*         Unit number on which to write the error at each iteration, */
/*         if this is desired for monitoring convergence.  If unit */
/*         number is 0, no writing will occur. */
/* R      :INOUT    Real R(N). */
/*         Work array used in calling routine.  It contains */
/*         information necessary to compute the residual RL = B-A*XL. */
/* Z      :WORK     Real Z(N). */
/*         Workspace used to hold the pseudo-residule M z = r. */
/* DZ     :WORK     Real DZ(N). */
/*         Workspace used to hold temporary vector(s). */
/* RWORK  :WORK     Real RWORK(USER DEFINABLE). */
/*         Real array that can be used by MSOLVE. */
/* IWORK  :WORK     Integer IWORK(USER DEFINABLE). */
/*         Integer array that can be used by MSOLVE. */
/* RNRM   :IN       Real. */
/*         Norm of the current residual.  Type of norm depends on ITOL. */
/* BNRM   :IN       Real. */
/*         Norm of the right hand side.  Type of norm depends on ITOL. */
/* SB     :IN       Real SB(N). */
/*         Scaling vector for B. */
/* SX     :IN       Real SX(N). */
/*         Scaling vector for X. */
/* JSCAL  :IN       Integer. */
/*         Flag indicating if scaling arrays SB and SX are being */
/*         used in the calling routine SPIGMR. */
/*         JSCAL=0 means SB and SX are not used and the */
/*                 algorithm will perform as if all */
/*                 SB(i) = 1 and SX(i) = 1. */
/*         JSCAL=1 means only SX is used, and the algorithm */
/*                 performs as if all SB(i) = 1. */
/*         JSCAL=2 means only SB is used, and the algorithm */
/*                 performs as if all SX(i) = 1. */
/*         JSCAL=3 means both SB and SX are used. */
/* KMP    :IN       Integer */
/*         The number of previous vectors the new vector VNEW */
/*         must be made orthogonal to.  (KMP .le. MAXL) */
/* LGMR   :IN       Integer */
/*         The number of GMRES iterations performed on the current call */
/*         to SPIGMR (i.e., # iterations since the last restart) and */
/*         the current order of the upper hessenberg */
/*         matrix HES. */
/* MAXL   :IN       Integer */
/*         The maximum allowable order of the matrix H. */
/* MAXLP1 :IN       Integer */
/*         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES. */
/* V      :IN       Real V(N,MAXLP1) */
/*         The N by (LGMR+1) array containing the LGMR */
/*         orthogonal vectors V(*,1) to V(*,LGMR). */
/* Q      :IN       Real Q(2*MAXL) */
/*         A real array of length 2*MAXL containing the components */
/*         of the Givens rotations used in the QR decomposition */
/*         of HES. */
/* SNORMW :IN       Real */
/*         A scalar containing the scaled norm of VNEW before it */
/*         is renormalized in SPIGMR. */
/* PROD   :IN       Real */
/*        The product s1*s2*...*sl = the product of the sines of the */
/*        givens rotations used in the QR factorization of */
/*        the hessenberg matrix HES. */
/* R0NRM  :IN       Real */
/*         The scaled norm of initial residual R0. */
/* HES    :IN       Real HES(MAXLP1,MAXL) */
/*         The upper triangular factor of the QR decomposition */
/*         of the (LGMR+1) by LGMR upper Hessenberg matrix whose */
/*         entries are the scaled inner-products of A*V(*,I) */
/*         and V(*,K). */
/* JPRE   :IN       Integer */
/*         Preconditioner type flag. */

/* *Description */
/*       When using the GMRES solver,  the preferred value  for ITOL */
/*       is 0.  This is due to the fact that when ITOL=0 the norm of */
/*       the residual required in the stopping test is  obtained for */
/*       free, since this value is already  calculated  in the GMRES */
/*       algorithm.   The  variable  RNRM contains the   appropriate */
/*       norm, which is equal to norm(SB*(RL - A*XL))  when right or */
/*       no   preconditioning is  being  performed,   and equal   to */
/*       norm(SB*Minv*(RL - A*XL))  when using left preconditioning. */
/*       Here, norm() is the Euclidean norm.  Nonzero values of ITOL */
/*       require  additional work  to  calculate the  actual  scaled */
/*       residual  or its scaled/preconditioned  form,  and/or   the */
/*       approximate solution XL.  Hence, these values of  ITOL will */
/*       not be as efficient as ITOL=0. */

/* ***ROUTINES CALLED     MSOLVE, SNRM2, SCOPY, */
/* ***END PROLOG ISSGMR */

/* ***FIRST EXECUTABLE STATEMENT ISSGMR */
    /* Parameter adjustments */
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --b;
    --x;
    --xl;
    --ia;
    --ja;
    --a;
    --r__;
    --z__;
    --dz;
    --rwork;
    --iwork;
    --sb;
    --sx;
    hes_dim1 = *maxlp1;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;
    --q;

    /* Function Body */
    ret_val = 0;
    if (*itol == 0) {

/*       Use input from SPIGMR to determine if stop conditions are met. */

	*err = *rnrm / *bnrm;
    }
    if (*itol > 0 && *itol <= 3) {

/*       Use SRLCAL to calculate the scaled residual vector. */
/*       Store answer in R. */

	if (*lgmr != 0) {
	    srlcal_(n, kmp, lgmr, maxl, &v[v_offset], &q[1], &r__[1], snormw, 
		    prod, r0nrm);
	}
	if (*itol <= 2) {
/*         err = ||Residual||/||RightHandSide||(2-Norms). */
	    *err = snrm2_(n, &r__[1], &c__1) / *bnrm;

/*         Unscale R by R0NRM*PROD when KMP < MAXL. */

	    if (*kmp < *maxl && *lgmr != 0) {
		tem = 1.f / (*r0nrm * *prod);
		sscal_(n, &tem, &r__[1], &c__1);
	    }
	} else if (*itol == 3) {
/*         err = Max |(Minv*Residual)(i)/x(i)| */
/*         When jpre .lt. 0, r already contains Minv*Residual. */
	    if (*jpre > 0) {
		(*msolve)(n, &r__[1], &dz[1], nelt, &ia[1], &ja[1], &a[1], 
			isym, &rwork[1], &iwork[1]);
		++(*nmsl);
	    }

/*         Unscale R by R0NRM*PROD when KMP < MAXL. */

	    if (*kmp < *maxl && *lgmr != 0) {
		tem = 1.f / (*r0nrm * *prod);
		sscal_(n, &tem, &r__[1], &c__1);
	    }

	    fuzz = r1mach_(&c__1);
	    ielmax = 1;
/* Computing MAX */
	    r__1 = dabs(x[1]);
	    ratmax = dabs(dz[1]) / dmax(r__1,fuzz);
	    i__1 = *n;
	    for (i__ = 2; i__ <= i__1; ++i__) {
/* Computing MAX */
		r__3 = (r__2 = x[i__], dabs(r__2));
		rat = (r__1 = dz[i__], dabs(r__1)) / dmax(r__3,fuzz);
		if (rat > ratmax) {
		    ielmax = i__;
		    ratmax = rat;
		}
/* L25: */
	    }
	    *err = ratmax;
	    if (ratmax <= *tol) {
		ret_val = 1;
	    }
	    if (*iunit > 0) {
		io___115.ciunit = *iunit;
		s_wsfe(&io___115);
		do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ielmax, (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&ratmax, (ftnlen)sizeof(real));
		e_wsfe();
	    }
	    return ret_val;
	}
    }
    if (*itol == 11) {

/*       Use SXLCAL to calculate the approximate solution XL. */

	if (*lgmr != 0 && *iter > 0) {
	    sxlcal_(n, lgmr, &x[1], &xl[1], &xl[1], &hes[hes_offset], maxlp1, 
		    &q[1], &v[v_offset], r0nrm, &dz[1], &sx[1], jscal, jpre, (
		    S_fp)msolve, nmsl, &rwork[1], &iwork[1], nelt, &ia[1], &
		    ja[1], &a[1], isym);
	} else if (*iter == 0) {
/*         Copy X to XL to check if initial guess is good enough. */
	    scopy_(n, &x[1], &c__1, &xl[1], &c__1);
	} else {
/*         Return since this is the first call to SPIGMR on a restart. */
	    return ret_val;
	}

	if (*jscal == 0 || *jscal == 2) {
/*         err = ||x-TrueSolution||/||TrueSolution||(2-Norms). */
	    if (*iter == 0) {
		solnrm = snrm2_(n, solblk_1.soln, &c__1);
	    }
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		dz[i__] = xl[i__] - solblk_1.soln[i__ - 1];
/* L30: */
	    }
	    *err = snrm2_(n, &dz[1], &c__1) / solnrm;
	} else {
	    if (*iter == 0) {
		solnrm = 0.f;
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		    r__1 = sx[i__] * solblk_1.soln[i__ - 1];
		    solnrm += r__1 * r__1;
/* L40: */
		}
		solnrm = sqrt(solnrm);
	    }
	    dxnrm = 0.f;
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		r__1 = sx[i__] * (xl[i__] - solblk_1.soln[i__ - 1]);
		dxnrm += r__1 * r__1;
/* L50: */
	    }
	    dxnrm = sqrt(dxnrm);
/*         err = ||SX*(x-TrueSolution)||/||SX*TrueSolution|| (2-Norms). */
	    *err = dxnrm / solnrm;
	}
    }

    if (*iunit != 0) {
	if (*iter == 0) {
	    io___118.ciunit = *iunit;
	    s_wsfe(&io___118);
	    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*itol), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*maxl), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*kmp), (ftnlen)sizeof(integer));
	    e_wsfe();
	}
	io___119.ciunit = *iunit;
	s_wsfe(&io___119);
	do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
	r__1 = *rnrm / *bnrm;
	do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
	do_fio(&c__1, (char *)&(*err), (ftnlen)sizeof(real));
	e_wsfe();
    }
    if (*err <= *tol) {
	ret_val = 1;
    }

    return ret_val;
/* ------------- LAST LINE OF ISSGMR FOLLOWS ---------------------------- */
} /* issgmr_ */

