/* slaputil.f -- translated by f2c (version 20100827).
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
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__0 = 0;
static real c_b50 = 0.f;
static integer c__59 = 59;
static integer c__55 = 55;

/* DECK SBHIN */
/* Subroutine */ int sbhin_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, real *soln, real *rhs, integer *iunit, 
	integer *job)
{
    /* Format strings */
    static char fmt_9000[] = "(a80)";
    static char fmt_9010[] = "(5i14)";
    static char fmt_9020[] = "(a3,11x,4i14)";
    static char fmt_9030[] = "(2a16,2a20)";

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void),
	     s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__, j;
    static char code[3];
    static integer ibgn, iend, nele, icol, nind, ncol;
    static real temp;
    static integer npls, nrow, nline, itemp;
    static char title[80];
    static integer nrils, nnvls, jobret;
    static char rinfmt[16], rhsfmt[20], nvlfmt[20], pntfmt[16];
    static integer nrhsls;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 0, 0, fmt_9000, 0 };
    static cilist io___3 = { 0, 0, 0, fmt_9010, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_9020, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_9030, 0 };
    static cilist io___21 = { 0, 0, 0, pntfmt, 0 };
    static cilist io___23 = { 0, 0, 0, rinfmt, 0 };
    static cilist io___24 = { 0, 0, 0, nvlfmt, 0 };
    static cilist io___25 = { 0, 5, 0, rhsfmt, 0 };


/* ***BEGIN PROLOGUE  SBHIN */
/* ***DATE WRITTEN   881107   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SBHIN-S), */
/*             Linear system, SLAP Sparse, Diagnostics */
/* ***AUTHOR  Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Read a Sparse Linear System in the Boeing/Harwell Format. */
/*            The matrix is read in and if the right hand side is also */
/*            present in the input file then it too is read in. */
/*            The matrix is then modified to be in the SLAP Column */
/*            format. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB */
/*     REAL    A(NELT), SOLN(N), RHS(N) */

/*     CALL SBHIN( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB ) */

/* *Arguments: */
/* N      :OUT      Integer */
/*         Order of the Matrix. */
/* NELT   :INOUT    Integer. */
/*         On input NELT is the maximum number of non-zeros that */
/*         can be stored in the IA, JA, A arrays. */
/*         On output NELT is the number of non-zeros stored in A. */
/* IA     :OUT      Integer IA(NELT). */
/* JA     :OUT      Integer JA(NELT). */
/* A      :OUT      Real A(NELT). */
/*         On output these arrays hold the matrix A in the SLAP */
/*         Triad format.  See "LONG DESCRIPTION", below. */
/* ISYM   :OUT      Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* SOLN   :OUT      Real SOLN(N). */
/*         The solution to the linear system, if present.  This array */
/*         is accessed if and only if JOB to read it in, see below. */
/*         If the user requests that SOLN be read in, but it is not in */
/*         the file, then it is simply zeroed out. */
/* RHS    :OUT      Real RHS(N). */
/*         The right hand side vector.  This array is accessed if and */
/*         only if JOB is set to read it in, see below. */
/*         If the user requests that RHS be read in, but it is not in */
/*         the file, then it is simply zeroed out. */
/* IUNIT  :IN       Integer. */
/*         Fortran logical I/O device unit number to write the matrix */
/*         to.  This unit must be connected in a system dependent fashion */
/*         to a file or the console or you will get a nasty message */
/*         from the Fortran I/O libraries. */
/* JOB    :INOUT    Integer. */
/*         Flag indicating what I/O operations to perform. */
/*         On input JOB indicates what Input operations to try to */
/*         perform. */
/*         JOB = 0 => Read only the matrix. */
/*             = 1 => Read matrix and RHS (if present). */
/*             = 2 => Read matrix and SOLN (if present). */
/*             = 3 => Read matrix, RHS and SOLN (if present). */
/*         On output JOB indicates what operations were actually */
/*         performed. */
/*               -3 => Unable to parse matrix "CODE" from input file */
/*                     to determine if only the lower triangle of matrix */
/*                     is stored. */
/*               -2 => Number of non-zeros (NELT) too large. */
/*               -1 => System size (N) too large. */
/*         JOB =  0 => Read in only the matrix. */
/*             =  1 => Read in the matrix and RHS. */
/*             =  2 => Read in the matrix and SOLN. */
/*             =  3 => Read in the matrix, RHS and SOLN. */
/*             = 10 => Read in only the matrix *STRUCTURE*, but no */
/*                     non-zero entries.  Hence, A(*) is not referenced */
/*                     and has the return values the same as the input. */
/*             = 11 => Read in the matrix *STRUCTURE* and RHS. */
/*             = 12 => Read in the matrix *STRUCTURE* and SOLN. */
/*             = 13 => Read in the matrix *STRUCTURE*, RHS and SOLN. */

/* *Precision:           Single Precision */
/* *Portability: */
/*         You must make sure that IUNIT is a valid Fortran logical */
/*         I/O device unit number and that the unit number has been */
/*         associated with a file or the console.  This is a system */
/*         dependent function. */

/* ***LONG DESCRIPTION */
/*       The format for the output is as follows.  On  the first line */
/*       are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT */
/*       and ISYM are described above.  IRHS is  a flag indicating if */
/*       the RHS was  written out (1 is  yes, 0 is  no).  ISOLN  is a */
/*       flag indicating if the SOLN was written out  (1 is yes, 0 is */
/*       no).  The format for the fist line is: 5i10.  Then comes the */
/*       NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format */
/*       for  these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then  comes */
/*       RHS(I), I = 1, N, if IRHS = 1.  Then  comes SOLN(I), I  = 1, */
/*       N, if ISOLN = 1.  The format for these lines is: 1X,E16.7. */

/*       =================== S L A P Triad format =================== */
/*       This routine requires that the  matrix A be   stored in  the */
/*       SLAP  Triad format.  In  this format only the non-zeros  are */
/*       stored.  They may appear in  *ANY* order.  The user supplies */
/*       three arrays of  length NELT, where  NELT is  the number  of */
/*       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For */
/*       each non-zero the user puts the row and column index of that */
/*       matrix element  in the IA and  JA arrays.  The  value of the */
/*       non-zero  matrix  element is  placed   in  the corresponding */
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
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SBHIN */

/*         Local Variables */



/*         Read Matrices In BOEING-HARWELL format. */

/* NLINE  Number of Data (after the header) lines in the file. */
/* NPLS   Number of lines for the Column Pointer data in the file. */
/* NRILS  Number of lines for the Row indicies in the data file. */
/* NNVLS  Number of lines for the Matrix elements in the data file. */
/* NRHSLS Number of lines for the RHS in the data file. */

/* ***FIRST EXECUTABLE STATEMENT  SBHIN */
    /* Parameter adjustments */
    --rhs;
    --soln;
    --a;
    --ja;
    --ia;

    /* Function Body */
    io___1.ciunit = *iunit;
    s_rsfe(&io___1);
    do_fio(&c__1, title, (ftnlen)80);
    e_rsfe();
    io___3.ciunit = *iunit;
    s_rsfe(&io___3);
    do_fio(&c__1, (char *)&nline, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&npls, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nrils, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nnvls, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nrhsls, (ftnlen)sizeof(integer));
    e_rsfe();
    io___9.ciunit = *iunit;
    s_rsfe(&io___9);
    do_fio(&c__1, code, (ftnlen)3);
    do_fio(&c__1, (char *)&nrow, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ncol, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nind, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nele, (ftnlen)sizeof(integer));
    e_rsfe();
    io___15.ciunit = *iunit;
    s_rsfe(&io___15);
    do_fio(&c__1, pntfmt, (ftnlen)16);
    do_fio(&c__1, rinfmt, (ftnlen)16);
    do_fio(&c__1, nvlfmt, (ftnlen)20);
    do_fio(&c__1, rhsfmt, (ftnlen)20);
    e_rsfe();

    if (nrow > *n) {
	*n = nrow;
	jobret = -1;
	goto L999;
    }
    if (nind > *nelt) {
	*nelt = nind;
	jobret = -2;
	goto L999;
    }

/*         Set the parameters. */

    *n = nrow;
    *nelt = nind;
    if (s_cmp(code, "RUA", (ftnlen)3, (ftnlen)3) == 0) {
	*isym = 0;
    } else if (s_cmp(code, "RSA", (ftnlen)3, (ftnlen)3) == 0) {
	*isym = 1;
    } else {
	jobret = -3;
	goto L999;
    }
    io___21.ciunit = *iunit;
    s_rsfe(&io___21);
    i__1 = *n + 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ja[i__], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    io___23.ciunit = *iunit;
    s_rsfe(&io___23);
    i__1 = *nelt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&ia[i__], (ftnlen)sizeof(integer));
    }
    e_rsfe();
    jobret = 10;
    if (nnvls > 0) {
	io___24.ciunit = *iunit;
	s_rsfe(&io___24);
	i__1 = *nelt;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&a[i__], (ftnlen)sizeof(real));
	}
	e_rsfe();
	jobret = 0;
    }
    if (nrhsls > 0 && *job % 2 == 1) {
	s_rsfe(&io___25);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&rhs[i__], (ftnlen)sizeof(real));
	}
	e_rsfe();
	++jobret;
    }

/*         Now loop thru the IA(i) array making sure that the Diagonal */
/*         matrix element appears first in the column.  Then sort the */
/*         rest of the column in ascending order. */

/* VD$R NOCONCUR */
/* VD$R NOVECTOR */
    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	ibgn = ja[icol];
	iend = ja[icol + 1] - 1;
	i__2 = iend;
	for (i__ = ibgn; i__ <= i__2; ++i__) {
	    if (ia[i__] == icol) {
/*         Swap the diag element with the first element in the column. */
		itemp = ia[i__];
		ia[i__] = ia[ibgn];
		ia[ibgn] = itemp;
		temp = a[i__];
		a[i__] = a[ibgn];
		a[ibgn] = temp;
		goto L40;
	    }
/* L30: */
	}
L40:
	++ibgn;
	if (ibgn < iend) {
	    i__2 = iend;
	    for (i__ = ibgn; i__ <= i__2; ++i__) {
		i__3 = iend;
		for (j = i__ + 1; j <= i__3; ++j) {
		    if (ia[i__] > ia[j]) {
			itemp = ia[i__];
			ia[i__] = ia[j];
			ia[j] = itemp;
			temp = a[i__];
			a[i__] = a[j];
			a[j] = temp;
		    }
/* L50: */
		}
/* L60: */
	    }
	}
/* L70: */
    }

/*         Set return flag. */
L999:
    *job = jobret;
    return 0;
/* ------------- LAST LINE OF SBHIN FOLLOWS ------------------------------ */
} /* sbhin_ */

/* DECK SCHKW */
/* Subroutine */ int schkw_(char *name__, integer *lociw, integer *leniw, 
	integer *locw, integer *lenw, integer *ierr, integer *iter, real *err,
	 ftnlen name_len)
{
    /* System generated locals */
    address a__1[3];
    integer i__1[3], i__2;

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    integer i_len(char *, ftnlen);

    /* Local variables */
    static char mesg[72];
    extern doublereal r1mach_(integer *);
    extern /* Subroutine */ int xerrwv_(char *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, real *, real *, 
	    ftnlen);

/* ***BEGIN PROLOGUE  SCHKW */
/* ***DATE WRITTEN   880225   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  R2 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SCHKW-S), */
/*             SLAP, Error Checking, Workspace Checking */
/* ***AUTHOR  Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP WORK/IWORK Array Bounds Checker. */
/*            This routine checks the work array lengths  and  inter- */
/*            faces to the SLATEC  error  handler  if  a  problem  is */
/*            found. */
/* ***DESCRIPTION */
/* *Usage: */
/*     CHARACTER*(*) NAME */
/*     INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER */
/*     REAL    ERR */

/*     CALL SCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR ) */

/* *Arguments: */
/* NAME   :IN       Character*(*). */
/*         Name of the calling routine.  This is used in the output */
/*         message, if an error is detected. */
/* LOCIW  :IN       Integer. */
/*         Location of the first free element in the integer workspace */
/*         array. */
/* LENIW  :IN       Integer. */
/*         Length of the integer workspace array. */
/* LOCW   :IN       Integer. */
/*         Location of the first free element in the real workspace */
/*         array. */
/* LENRW  :IN       Integer. */
/*         Length of the real workspace array. */
/* IERR   :OUT      Integer. */
/*         Return error flag. */
/*               IERR = 0 => All went well. */
/*               IERR = 1 => Insufficient storage allocated for */
/*                           WORK or IWORK. */
/* ITER   :OUT      Integer. */
/*         Set to 0 if an error is detected. */
/* ERR    :OUT      Real. */
/*         Set to a very large number if an error is detected. */

/* *Precision:           Single Precision */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  R1MACH, XERRWV */
/* ***END PROLOGUE  SCHKW */

/*         Check the Integer workspace situation. */
/* ***FIRST EXECUTABLE STATEMENT  SCHKW */
    *ierr = 0;
    if (*lociw > *leniw) {
	*ierr = 1;
	*iter = 0;
	*err = r1mach_(&c__2);
/* Writing concatenation */
	i__1[0] = name_len, a__1[0] = name__;
	i__1[1] = 32, a__1[1] = ": INTEGER work array too short. ";
	i__1[2] = 35, a__1[2] = " IWORK needs i1: have allocated i2.";
	s_cat(mesg, a__1, i__1, &c__3, (ftnlen)72);
	i__2 = i_len(mesg, (ftnlen)72);
	xerrwv_(mesg, &i__2, &c__1, &c__1, &c__2, lociw, leniw, &c__0, &c_b50,
		 &c_b50, (ftnlen)72);
    }

/*         Check the Real workspace situation. */
    if (*locw > *lenw) {
	*ierr = 1;
	*iter = 0;
	*err = r1mach_(&c__2);
/* Writing concatenation */
	i__1[0] = name_len, a__1[0] = name__;
	i__1[1] = 29, a__1[1] = ": REAL work array too short. ";
	i__1[2] = 35, a__1[2] = " RWORK needs i1: have allocated i2.";
	s_cat(mesg, a__1, i__1, &c__3, (ftnlen)72);
	i__2 = i_len(mesg, (ftnlen)72);
	xerrwv_(mesg, &i__2, &c__1, &c__1, &c__2, locw, lenw, &c__0, &c_b50, &
		c_b50, (ftnlen)72);
    }
    return 0;
/* ------------- LAST LINE OF SCHKW FOLLOWS ---------------------------- */
} /* schkw_ */

/* DECK QS2I1R */
/* Subroutine */ int qs2i1r_(integer *ia, integer *ja, real *a, integer *n, 
	integer *kflag)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, j, k, l, m;
    static real r__;
    static integer ij, il[21];
    static real ta;
    static integer kk, nn, it, iu[21], jt, iit, jjt;
    static real tta;
    extern /* Subroutine */ int xerror_(char *, integer *, integer *, integer 
	    *, ftnlen);

/* ***BEGIN PROLOGUE  QS2I1R */
/* ***DATE WRITTEN   761118   (YYMMDD) */
/* ***REVISION DATE  890125   (YYMMDD) */
/* ***CATEGORY NO.  N6A2A */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=INTEGER(QS2I1R-I), */
/*             QUICKSORT,SINGLETON QUICKSORT,SORT,SORTING */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/*           Kahaner, D. K., (NBS) */
/*           Seager, M. K., (LLNL) seager@lll-crg.llnl.gov */
/*           Wisniewski, J. A., (SNLA) */
/* ***PURPOSE  Sort an integer array also moving an integer and real array */
/*            This routine sorts the integer array IA and makes the same */
/*            interchanges in the integer array JA and the real array A. */
/*            The array IA may be sorted in increasing order or decreas- */
/*            ing order. A slightly modified QUICKSORT algorithm is used. */

/* ***DESCRIPTION */
/*     Written by Rondall E Jones */
/*     Modified by John A. Wisniewski to use the Singleton QUICKSORT */
/*     algorithm. date 18 November 1976. */

/*     Further modified by David K. Kahaner */
/*     National Bureau of Standards */
/*     August, 1981 */

/*     Even further modification made to bring the code up to the */
/*     Fortran 77 level and make it more readable and to carry */
/*     along one integer array and one real array during the sort by */
/*     Mark K. Seager */
/*     Lawrence Livermore National Laboratory */
/*     November, 1987 */
/*     This routine was adapted from the ISORT routine. */

/*     ABSTRACT */
/*         This routine sorts an integer array IA and makes the same */
/*         interchanges in the integer array JA and the real array A. */
/*         The array a may be sorted in increasing order or decreasing */
/*         order.  A slightly modified quicksort algorithm is used. */

/*     DESCRIPTION OF PARAMETERS */
/*        IA - Integer array of values to be sorted. */
/*        JA - Integer array to be carried along. */
/*         A - Real array to be carried along. */
/*         N - Number of values in integer array IA to be sorted. */
/*     KFLAG - Control parameter */
/*           = 1 means sort IA in INCREASING order. */
/*           =-1 means sort IA in DECREASING order. */

/* ***REFERENCES */
/*     Singleton, R. C., Algorithm 347, "An Efficient Algorithm for */
/*     Sorting with Minimal Storage", cacm, Vol. 12, No. 3, 1969, */
/*     Pp. 185-187. */
/* ***ROUTINES CALLED  XERROR */
/* ***END PROLOGUE  QS2I1R */
/* VD$R NOVECTOR */
/* VD$R NOCONCUR */

/* ***FIRST EXECUTABLE STATEMENT  QS2I1R */
    /* Parameter adjustments */
    --a;
    --ja;
    --ia;

    /* Function Body */
    nn = *n;
    if (nn < 1) {
	xerror_("QS2I1R- the number of values to be sorted was not POSITIVE.",
		 &c__59, &c__1, &c__1, (ftnlen)59);
	return 0;
    }
    if (*n == 1) {
	return 0;
    }
    kk = abs(*kflag);
    if (kk != 1) {
	xerror_("QS2I1R- the sort control parameter, k, was not 1 OR -1.", &
		c__55, &c__2, &c__1, (ftnlen)55);
	return 0;
    }

/*     Alter array IA to get decreasing order if needed. */

    if (*kflag < 1) {
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia[i__] = -ia[i__];
/* L20: */
	}
    }

/*     Sort IA and carry JA and A along. */
/*     And now...Just a little black magic... */
    m = 1;
    i__ = 1;
    j = nn;
    r__ = .375f;
L210:
    if (r__ <= .5898437f) {
	r__ += .0390625f;
    } else {
	r__ += -.21875f;
    }
L225:
    k = i__;

/*     Select a central element of the array and save it in location */
/*     it, jt, at. */

    ij = i__ + (integer) ((real) (j - i__) * r__);
    it = ia[ij];
    jt = ja[ij];
    ta = a[ij];

/*     If first element of array is greater than it, interchange with it. */

    if (ia[i__] > it) {
	ia[ij] = ia[i__];
	ia[i__] = it;
	it = ia[ij];
	ja[ij] = ja[i__];
	ja[i__] = jt;
	jt = ja[ij];
	a[ij] = a[i__];
	a[i__] = ta;
	ta = a[ij];
    }
    l = j;

/*     If last element of array is less than it, swap with it. */

    if (ia[j] < it) {
	ia[ij] = ia[j];
	ia[j] = it;
	it = ia[ij];
	ja[ij] = ja[j];
	ja[j] = jt;
	jt = ja[ij];
	a[ij] = a[j];
	a[j] = ta;
	ta = a[ij];

/*     If first element of array is greater than it, swap with it. */

	if (ia[i__] > it) {
	    ia[ij] = ia[i__];
	    ia[i__] = it;
	    it = ia[ij];
	    ja[ij] = ja[i__];
	    ja[i__] = jt;
	    jt = ja[ij];
	    a[ij] = a[i__];
	    a[i__] = ta;
	    ta = a[ij];
	}
    }

/*     Find an element in the second half of the array which is */
/*     smaller than it. */

L240:
    --l;
    if (ia[l] > it) {
	goto L240;
    }

/*     Find an element in the first half of the array which is */
/*     greater than it. */

L245:
    ++k;
    if (ia[k] < it) {
	goto L245;
    }

/*     Interchange these elements. */

    if (k <= l) {
	iit = ia[l];
	ia[l] = ia[k];
	ia[k] = iit;
	jjt = ja[l];
	ja[l] = ja[k];
	ja[k] = jjt;
	tta = a[l];
	a[l] = a[k];
	a[k] = tta;
	goto L240;
    }

/*     Save upper and lower subscripts of the array yet to be sorted. */

    if (l - i__ > j - k) {
	il[m - 1] = i__;
	iu[m - 1] = l;
	i__ = k;
	++m;
    } else {
	il[m - 1] = k;
	iu[m - 1] = j;
	j = l;
	++m;
    }
    goto L260;

/*     Begin again on another portion of the unsorted array. */

L255:
    --m;
    if (m == 0) {
	goto L300;
    }
    i__ = il[m - 1];
    j = iu[m - 1];
L260:
    if (j - i__ >= 1) {
	goto L225;
    }
    if (i__ == j) {
	goto L255;
    }
    if (i__ == 1) {
	goto L210;
    }
    --i__;
L265:
    ++i__;
    if (i__ == j) {
	goto L255;
    }
    it = ia[i__ + 1];
    jt = ja[i__ + 1];
    ta = a[i__ + 1];
    if (ia[i__] <= it) {
	goto L265;
    }
    k = i__;
L270:
    ia[k + 1] = ia[k];
    ja[k + 1] = ja[k];
    a[k + 1] = a[k];
    --k;
    if (it < ia[k]) {
	goto L270;
    }
    ia[k + 1] = it;
    ja[k + 1] = jt;
    a[k + 1] = ta;
    goto L265;

/*     Clean up, if necessary. */

L300:
    if (*kflag < 1) {
	i__1 = nn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    ia[i__] = -ia[i__];
/* L310: */
	}
    }
    return 0;
/* ------------- LAST LINE OF QS2I1R FOLLOWS ---------------------------- */
} /* qs2i1r_ */

/* DECK SS2Y */
/* Subroutine */ int ss2y_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, ibgn, iend, icol;
    static real temp;
    static integer itemp;
    extern /* Subroutine */ int qs2i1r_(integer *, integer *, real *, integer 
	    *, integer *);

/* ***BEGIN PROLOGUE  SS2Y */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SS2Y-S), */
/*             Linear system, SLAP Sparse */
/* ***AUTHOR  Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  SLAP Triad to SLAP Column Format Converter. */
/*            Routine to convert from the SLAP Triad to SLAP Column */
/*            format. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM */
/*     REAL    A(NELT) */

/*     CALL SS2Y( N, NELT, IA, JA, A, ISYM ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of non-zeros stored in A. */
/* IA     :INOUT    Integer IA(NELT). */
/* JA     :INOUT    Integer JA(NELT). */
/* A      :INOUT    Real A(NELT). */
/*         These arrays should hold the matrix A in either the SLAP */
/*         Triad format or the SLAP Column format.  See "LONG */
/*         DESCRIPTION", below.  If the SLAP Triad format is used */
/*         this format is translated to the SLAP Column format by */
/*         this routine. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */

/* *Precision:           Single Precision */

/* ***LONG DESCRIPTION */
/*       The Sparse Linear Algebra Package (SLAP) utilizes two matrix */
/*       data structures: 1) the  SLAP Triad  format or  2)  the SLAP */
/*       Column format.  The user can hand this routine either of the */
/*       of these data structures.  If the SLAP Triad format is give */
/*       as input then this routine transforms it into SLAP Column */
/*       format.  The way this routine tells which format is given as */
/*       input is to look at JA(N+1).  If JA(N+1) = NELT+1 then we */
/*       have the SLAP Column format.  If that equality does not hold */
/*       then it is assumed that the IA, JA, A arrays contain the */
/*       SLAP Triad format. */

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

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  QS2I1R */
/* ***END PROLOGUE  SS2Y */

/*         Check to see if the (IA,JA,A) arrays are in SLAP Column */
/*         format.  If it's not then transform from SLAP Triad. */
/* ***FIRST EXECUTABLE STATEMENT  SS2LT */
    /* Parameter adjustments */
    --a;
    --ja;
    --ia;

    /* Function Body */
    if (ja[*n + 1] == *nelt + 1) {
	return 0;
    }

/*         Sort into ascending order by COLUMN (on the ja array). */
/*         This will line up the columns. */

    qs2i1r_(&ja[1], &ia[1], &a[1], nelt, &c__1);

/*         Loop over each column to see where the column indicies change */
/*         in the column index array ja.  This marks the beginning of the */
/*         next column. */

/* VD$R NOVECTOR */
    ja[1] = 1;
    i__1 = *n - 1;
    for (icol = 1; icol <= i__1; ++icol) {
	i__2 = *nelt;
	for (j = ja[icol] + 1; j <= i__2; ++j) {
	    if (ja[j] != icol) {
		ja[icol + 1] = j;
		goto L20;
	    }
/* L10: */
	}
L20:
	;
    }
    ja[*n + 1] = *nelt + 1;

/*         Mark the n+2 element so that future calls to a SLAP routine */
/*         utilizing the YSMP-Column storage format will be able to tell. */

    ja[*n + 2] = 0;

/*         Now loop thru the ia(i) array making sure that the Diagonal */
/*         matrix element appears first in the column.  Then sort the */
/*         rest of the column in ascending order. */

    i__1 = *n;
    for (icol = 1; icol <= i__1; ++icol) {
	ibgn = ja[icol];
	iend = ja[icol + 1] - 1;
	i__2 = iend;
	for (i__ = ibgn; i__ <= i__2; ++i__) {
	    if (ia[i__] == icol) {
/*         Swap the diag element with the first element in the column. */
		itemp = ia[i__];
		ia[i__] = ia[ibgn];
		ia[ibgn] = itemp;
		temp = a[i__];
		a[i__] = a[ibgn];
		a[ibgn] = temp;
		goto L40;
	    }
/* L30: */
	}
L40:
	++ibgn;
	if (ibgn < iend) {
	    i__2 = iend;
	    for (i__ = ibgn; i__ <= i__2; ++i__) {
		i__3 = iend;
		for (j = i__ + 1; j <= i__3; ++j) {
		    if (ia[i__] > ia[j]) {
			itemp = ia[i__];
			ia[i__] = ia[j];
			ia[j] = itemp;
			temp = a[i__];
			a[i__] = a[j];
			a[j] = temp;
		    }
/* L50: */
		}
/* L60: */
	    }
	}
/* L70: */
    }
    return 0;
/* ------------- LAST LINE OF SS2Y FOLLOWS ---------------------------- */
} /* ss2y_ */

/* DECK SCPPLT */
/* Subroutine */ int scpplt_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, integer *iunit)
{
    /* Format strings */
    static char fmt_1000[] = "(/\002**** Picture of Column SLAP matrix follo"
	    "ws ****\002/\002 N, NELT and Density = \002,2i10,e16.7)";
    static char fmt_1010[] = "(4x,255(i1))";
    static char fmt_1020[] = "(1x,i3,a)";

    /* System generated locals */
    integer i__1, i__2;
    real r__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, jbgn, jend, icol, nmax, irow;
    static char chmat[225*225];

    /* Fortran I/O blocks */
    static cilist io___65 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___66 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_1020, 0 };


/* ***BEGIN PROLOGUE  SCPPLT */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(SCPPLT-S), */
/*             Linear system, SLAP Sparse, Diagnostics */
/* ***AUTHOR  Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Printer Plot of SLAP Column Format Matrix. */
/*            Routine to print out a SLAP Column format matrix in */
/*            a "printer plot" graphical representation. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(N+1), ISYM, IUNIT */
/*     REAL    A(NELT) */

/*     CALL SCPPLT( N, NELT, IA, JA, A, ISYM, IUNIT ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of non-zeros stored in A. */
/* IA     :INOUT    Integer IA(NELT). */
/* JA     :INOUT    Integer JA(N+1). */
/* A      :INOUT    Real A(NELT). */
/*         These arrays should hold the matrix A in the SLAP */
/*         Column format.  See "LONG DESCRIPTION", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* IUNIT  :IN       Integer. */
/*         Fortran logical I/O device unit number to write the matrix */
/*         to.  This unit must be connected in a system dependent fashion */
/*         to a file or the console or you will get a nasty message */
/*         from the Fortran I/O libraries. */

/* *Precision:           Single Precision */
/* *Portability: */
/*         You must make sure that IUNIT is a valid Fortran logical */
/*         I/O device unit number and that the unit number has been */
/*         associated with a file or the console.  This is a system */
/*         dependent function. */

/* ***LONG DESCRIPTION */
/*       This routine prints out a SLAP  Column format matrix  to the */
/*       Fortran logical I/O unit   number  IUNIT.  The  numbers them */
/*       selves  are not printed  out, but   rather  a one  character */
/*       representation of the numbers.   Elements of the matrix that */
/*       are not represented in the (IA,JA,A)  arrays are  denoted by */
/*       ' ' character (a blank).  Elements of A that are *ZERO* (and */
/*       hence  should  really not be  stored) are  denoted  by a '0' */
/*       character.  Elements of A that are *POSITIVE* are denoted by */
/*       'D' if they are Diagonal elements  and '#' if  they are off */
/*       Diagonal  elements.  Elements of  A that are *NEGATIVE* are */
/*       denoted by 'N'  if they  are Diagonal  elements and  '*' if */
/*       they are off Diagonal elements. */

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

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  SCPPLT */

/*         Set up the character matrix... */
/* ***FIRST EXECUTABLE STATEMENT  SCPPLT */
    /* Parameter adjustments */
    --a;
    --ja;
    --ia;

    /* Function Body */
    nmax = min(225,*n);
    i__1 = nmax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	s_copy(chmat + (i__ - 1) * 225, " ", nmax, (ftnlen)1);
/* L10: */
    }
    i__1 = nmax;
    for (icol = 1; icol <= i__1; ++icol) {
	jbgn = ja[icol];
	jend = ja[icol + 1] - 1;
	i__2 = jend;
	for (j = jbgn; j <= i__2; ++j) {
	    irow = ia[j];
	    if (irow <= nmax) {
		if (*isym != 0) {
/*         Put in non-sym part as well... */
		    if (a[j] == 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '0';
		    } else if (a[j] > 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '#';
		    } else {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '*';
		    }
		}
		if (irow == icol) {
/*         Diagonal entry. */
		    if (a[j] == 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '0';
		    } else if (a[j] > 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = 'D';
		    } else {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = 'N';
		    }
		} else {
/*         Off-Diagonal entry */
		    if (a[j] == 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '0';
		    } else if (a[j] > 0.f) {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '#';
		    } else {
			*(unsigned char *)&chmat[(irow - 1) * 225 + (icol - 1)
				] = '*';
		    }
		}
	    }
/* L20: */
	}
/* L30: */
    }

/*         Write out the heading. */
    io___65.ciunit = *iunit;
    s_wsfe(&io___65);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*nelt), (ftnlen)sizeof(integer));
    r__1 = (real) (*nelt) / (real) (*n * *n);
    do_fio(&c__1, (char *)&r__1, (ftnlen)sizeof(real));
    e_wsfe();
    io___66.ciunit = *iunit;
    s_wsfe(&io___66);
    i__1 = nmax;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ % 10;
	do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
    }
    e_wsfe();

/*         Write out the character representations matrix elements. */
    i__2 = nmax;
    for (irow = 1; irow <= i__2; ++irow) {
	io___67.ciunit = *iunit;
	s_wsfe(&io___67);
	do_fio(&c__1, (char *)&irow, (ftnlen)sizeof(integer));
	do_fio(&c__1, chmat + (irow - 1) * 225, nmax);
	e_wsfe();
/* L40: */
    }
    return 0;
/* ------------- LAST LINE OF SCPPLT FOLLOWS ---------------------------- */
} /* scpplt_ */

/* DECK STOUT */
/* Subroutine */ int stout_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, real *soln, real *rhs, integer *iunit, 
	integer *job)
{
    /* Format strings */
    static char fmt_1000[] = "(5i10)";
    static char fmt_1010[] = "(1x,i5,1x,i5,1x,e16.7)";
    static char fmt_1020[] = "(1x,e16.7)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, irhs, isoln;

    /* Fortran I/O blocks */
    static cilist io___70 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___72 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___73 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___74 = { 0, 0, 0, fmt_1020, 0 };


/* ***BEGIN PROLOGUE  STOUT */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(STOUT-S), */
/*             Linear system, SLAP Sparse, Diagnostics */
/* ***AUTHOR  Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Write out SLAP Triad Format Linear System. */
/*            Routine to write out a SLAP Triad format matrix and */
/*            right hand side and solution to the system, if known. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB */
/*     REAL    A(NELT), SOLN(N), RHS(N) */

/*     CALL STOUT( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB ) */

/* *Arguments: */
/* N      :IN       Integer */
/*         Order of the Matrix. */
/* NELT   :IN       Integer. */
/*         Number of non-zeros stored in A. */
/* IA     :INOUT    Integer IA(NELT). */
/* JA     :INOUT    Integer JA(NELT). */
/* A      :INOUT    Real A(NELT). */
/*         These arrays should hold the matrix A in the SLAP */
/*         Triad format.  See "LONG DESCRIPTION", below. */
/* ISYM   :IN       Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* SOLN   :IN       Real SOLN(N). */
/*         The solution to the linear system, if known.  This array */
/*         is accessed if and only if JOB is set to print it out, */
/*         see below. */
/* RHS    :IN       Real RHS(N). */
/*         The right hand side vector.  This array is accessed if and */
/*         only if JOB is set to print it out, see below. */
/* IUNIT  :IN       Integer. */
/*         Fortran logical I/O device unit number to write the matrix */
/*         to.  This unit must be connected in a system dependent fashion */
/*         to a file or the console or you will get a nasty message */
/*         from the Fortran I/O libraries. */
/* JOB    :IN       Integer. */
/*         Flag indicating what I/O operations to perform. */
/*         JOB = 0 => Print only the matrix. */
/*             = 1 => Print matrix and RHS. */
/*             = 2 => Print matrix and SOLN. */
/*             = 3 => Print matrix, RHS and SOLN. */

/* *Precision:           Single Precision */
/* *Portability: */
/*         You must make sure that IUNIT is a valid Fortran logical */
/*         I/O device unit number and that the unit number has been */
/*         associated with a file or the console.  This is a system */
/*         dependent function. */

/* ***LONG DESCRIPTION */
/*       The format for the output is as follows.  On  the first line */
/*       are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT */
/*       and ISYM are described above.  IRHS is  a flag indicating if */
/*       the RHS was  written out (1 is  yes, 0 is  no).  ISOLN  is a */
/*       flag indicating if the SOLN was written out  (1 is yes, 0 is */
/*       no).  The format for the fist line is: 5i10.  Then comes the */
/*       NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format */
/*       for  these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then  comes */
/*       RHS(I), I = 1, N, if IRHS = 1.  Then  comes SOLN(I), I  = 1, */
/*       N, if ISOLN = 1.  The format for these lines is: 1X,E16.7. */

/*       =================== S L A P Triad format =================== */
/*       This routine requires that the  matrix A be   stored in  the */
/*       SLAP  Triad format.  In  this format only the non-zeros  are */
/*       stored.  They may appear in  *ANY* order.  The user supplies */
/*       three arrays of  length NELT, where  NELT is  the number  of */
/*       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For */
/*       each non-zero the user puts the row and column index of that */
/*       matrix element  in the IA and  JA arrays.  The  value of the */
/*       non-zero  matrix  element is  placed   in  the corresponding */
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
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  STOUT */

/*         Local variables. */


/*         If RHS and SOLN are to be printed also. */
/*         Write out the information heading. */
/* ***FIRST EXECUTABLE STATEMENT  STOUT */
    /* Parameter adjustments */
    --rhs;
    --soln;
    --a;
    --ja;
    --ia;

    /* Function Body */
    irhs = 0;
    isoln = 0;
    if (*job == 1 || *job == 3) {
	irhs = 1;
    }
    if (*job > 1) {
	isoln = 1;
    }
    io___70.ciunit = *iunit;
    s_wsfe(&io___70);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*nelt), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*isym), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&irhs, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&isoln, (ftnlen)sizeof(integer));
    e_wsfe();

/*         Write out the matrix non-zeros in Triad format. */
    i__1 = *nelt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___72.ciunit = *iunit;
	s_wsfe(&io___72);
	do_fio(&c__1, (char *)&ia[i__], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ja[i__], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&a[i__], (ftnlen)sizeof(real));
	e_wsfe();
/* L10: */
    }

/*         If requested, write out the rhs. */
    if (irhs == 1) {
	io___73.ciunit = *iunit;
	s_wsfe(&io___73);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&rhs[i__], (ftnlen)sizeof(real));
	}
	e_wsfe();
    }

/*         If requested, write out the soln. */
    if (isoln == 1) {
	io___74.ciunit = *iunit;
	s_wsfe(&io___74);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&soln[i__], (ftnlen)sizeof(real));
	}
	e_wsfe();
    }
    return 0;
/* ------------- LAST LINE OF STOUT FOLLOWS ---------------------------- */
} /* stout_ */

/* DECK STIN */
/* Subroutine */ int stin_(integer *n, integer *nelt, integer *ia, integer *
	ja, real *a, integer *isym, real *soln, real *rhs, integer *iunit, 
	integer *job)
{
    /* Format strings */
    static char fmt_1000[] = "(5i10)";
    static char fmt_1010[] = "(1x,i5,1x,i5,1x,e16.7)";
    static char fmt_1020[] = "(1x,e16.7)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(void);

    /* Local variables */
    static integer i__, irhs, isoln, jobret, neltmax;

    /* Fortran I/O blocks */
    static cilist io___76 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___80 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___82 = { 0, 0, 0, fmt_1020, 0 };
    static cilist io___83 = { 0, 0, 0, fmt_1020, 0 };


/* ***BEGIN PROLOGUE  STIN */
/* ***DATE WRITTEN   871119   (YYMMDD) */
/* ***REVISION DATE  881213   (YYMMDD) */
/* ***CATEGORY NO.  D2A4, D2B4 */
/* ***KEYWORDS  LIBRARY=SLATEC(SLAP), */
/*             TYPE=SINGLE PRECISION(STIN-S), */
/*             Linear system, SLAP Sparse, Diagnostics */
/* ***AUTHOR  Seager, Mark K., (LLNL) */
/*             Lawrence Livermore National Laboratory */
/*             PO BOX 808, L-300 */
/*             Livermore, CA 94550 (415) 423-3141 */
/*             seager@lll-crg.llnl.gov */
/* ***PURPOSE  Read in SLAP Triad Format Linear System. */
/*            Routine to read in a SLAP Triad format matrix and */
/*            right hand side and solution to the system, if known. */
/* ***DESCRIPTION */
/* *Usage: */
/*     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IUNIT, JOB */
/*     REAL    A(NELT), SOLN(N), RHS(N) */

/*     CALL STIN( N, NELT, IA, JA, A, ISYM, SOLN, RHS, IUNIT, JOB ) */

/* *Arguments: */
/* N      :OUT      Integer */
/*         Order of the Matrix. */
/* NELT   :INOUT    Integer. */
/*         On input NELT is the maximum number of non-zeros that */
/*         can be stored in the IA, JA, A arrays. */
/*         On output NELT is the number of non-zeros stored in A. */
/* IA     :OUT      Integer IA(NELT). */
/* JA     :OUT      Integer JA(NELT). */
/* A      :OUT      Real A(NELT). */
/*         On output these arrays hold the matrix A in the SLAP */
/*         Triad format.  See "LONG DESCRIPTION", below. */
/* ISYM   :OUT      Integer. */
/*         Flag to indicate symmetric storage format. */
/*         If ISYM=0, all nonzero entries of the matrix are stored. */
/*         If ISYM=1, the matrix is symmetric, and only the lower */
/*         triangle of the matrix is stored. */
/* SOLN   :OUT      Real SOLN(N). */
/*         The solution to the linear system, if present.  This array */
/*         is accessed if and only if JOB to read it in, see below. */
/*         If the user requests that SOLN be read in, but it is not in */
/*         the file, then it is simply zeroed out. */
/* RHS    :OUT      Real RHS(N). */
/*         The right hand side vector.  This array is accessed if and */
/*         only if JOB is set to read it in, see below. */
/*         If the user requests that RHS be read in, but it is not in */
/*         the file, then it is simply zeroed out. */
/* IUNIT  :IN       Integer. */
/*         Fortran logical I/O device unit number to write the matrix */
/*         to.  This unit must be connected in a system dependent fashion */
/*         to a file or the console or you will get a nasty message */
/*         from the Fortran I/O libraries. */
/* JOB    :INOUT    Integer. */
/*         Flag indicating what I/O operations to perform. */
/*         On input JOB indicates what Input operations to try to */
/*         perform. */
/*         JOB = 0 => Read only the matrix. */
/*             = 1 => Read matrix and RHS (if present). */
/*             = 2 => Read matrix and SOLN (if present). */
/*             = 3 => Read matrix, RHS and SOLN (if present). */
/*         On output JOB indicates what operations were actually */
/*         performed. */
/*         JOB = 0 => Read in only the matrix. */
/*             = 1 => Read in the matrix and RHS. */
/*             = 2 => Read in the matrix and SOLN. */
/*             = 3 => Read in the matrix, RHS and SOLN. */

/* *Precision:           Single Precision */
/* *Portability: */
/*         You must make sure that IUNIT is a valid Fortran logical */
/*         I/O device unit number and that the unit number has been */
/*         associated with a file or the console.  This is a system */
/*         dependent function. */

/* ***LONG DESCRIPTION */
/*       The format for the output is as follows.  On  the first line */
/*       are counters and flags: N, NELT, ISYM, IRHS, ISOLN.  N, NELT */
/*       and ISYM are described above.  IRHS is  a flag indicating if */
/*       the RHS was  written out (1 is  yes, 0 is  no).  ISOLN  is a */
/*       flag indicating if the SOLN was written out  (1 is yes, 0 is */
/*       no).  The format for the fist line is: 5i10.  Then comes the */
/*       NELT Triad's IA(I), JA(I) and A(I), I = 1, NELT.  The format */
/*       for  these lines is   :  1X,I5,1X,I5,1X,E16.7.   Then  comes */
/*       RHS(I), I = 1, N, if IRHS = 1.  Then  comes SOLN(I), I  = 1, */
/*       N, if ISOLN = 1.  The format for these lines is: 1X,E16.7. */

/*       =================== S L A P Triad format =================== */
/*       This routine requires that the  matrix A be   stored in  the */
/*       SLAP  Triad format.  In  this format only the non-zeros  are */
/*       stored.  They may appear in  *ANY* order.  The user supplies */
/*       three arrays of  length NELT, where  NELT is  the number  of */
/*       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For */
/*       each non-zero the user puts the row and column index of that */
/*       matrix element  in the IA and  JA arrays.  The  value of the */
/*       non-zero  matrix  element is  placed   in  the corresponding */
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
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***END PROLOGUE  STIN */

/*         Local variables. */


/*         Read in the information heading. */
/* ***FIRST EXECUTABLE STATEMENT  STIN */
    /* Parameter adjustments */
    --rhs;
    --soln;
    --a;
    --ja;
    --ia;

    /* Function Body */
    neltmax = *nelt;
    io___76.ciunit = *iunit;
    s_rsfe(&io___76);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*nelt), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*isym), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&irhs, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&isoln, (ftnlen)sizeof(integer));
    e_rsfe();
    *nelt = min(*nelt,neltmax);

/*         Read in the matrix non-zeros in Triad format. */
    i__1 = *nelt;
    for (i__ = 1; i__ <= i__1; ++i__) {
	io___80.ciunit = *iunit;
	s_rsfe(&io___80);
	do_fio(&c__1, (char *)&ia[i__], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ja[i__], (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&a[i__], (ftnlen)sizeof(real));
	e_rsfe();
/* L10: */
    }

/*         If requested, read in the rhs. */
    jobret = 0;
    if (*job == 1 || *job == 3) {

/*         Check to see if rhs is in the file. */
	if (irhs == 1) {
	    jobret = 1;
	    io___82.ciunit = *iunit;
	    s_rsfe(&io___82);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&rhs[i__], (ftnlen)sizeof(real));
	    }
	    e_rsfe();
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		rhs[i__] = 0.f;
/* L20: */
	    }
	}
    }

/*         If requested, read in the soln. */
    if (*job > 1) {

/*         Check to see if soln is in the file. */
	if (isoln == 1) {
	    jobret += 2;
	    io___83.ciunit = *iunit;
	    s_rsfe(&io___83);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&soln[i__], (ftnlen)sizeof(real));
	    }
	    e_rsfe();
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		soln[i__] = 0.f;
/* L30: */
	    }
	}
    }

    *job = jobret;
    return 0;
/* ------------- LAST LINE OF STIN FOLLOWS ---------------------------- */
} /* stin_ */

