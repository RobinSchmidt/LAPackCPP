*DECK SCGN
      SUBROUTINE SCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, 
     $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, 
     $     ATP, ATZ, DZ, ATDZ, RWORK, IWORK)
C***BEGIN PROLOGUE  SCGN
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SCGN-S),
C             Non-Symmetric Linear system solve, Sparse, 
C             Iterative Precondition, Normal Equations.
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Preconditioned CG Sparse Ax=b Solver for Normal Equations.
C            Routine  to solve a general linear system Ax = b using the
C            Preconditioned Conjugate Gradient method  applied to   the
C            normal equations AA'y = b, x=A'y.
C***DESCRIPTION
C *Usage:
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C     INTEGER  ITER, IERR, IUNIT, IWORK(USER DEFINABLE)
C     REAL     B(N), X(N), A(NELT), TOL, ERR, R(N), Z(N), P(N), ATP(N)
C     REAL     ATZ(N), DZ(N), ATDZ(N)
C     REAL     RWORK(USER DEFINABLE)
C     EXTERNAL MATVEC, MTTVEC, MSOLVE
C
C     CALL SCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, 
C    $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, 
C    $     Z, P, ATP, ATZ, DZ, ATDZ, RWORK, IWORK)
C         
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right-hand side vector.
C X      :INOUT    Real X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Integer A(NELT).
C         These arrays contain the matrix data structure for A.
C         It could take any form.  See "Description", below
C         for more late breaking details...
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MATVEC :EXT      External.
C         Name of a routine which performs the matrix vector multiply
C         y = A*X given A and X.  The name of the MATVEC routine must 
C         be declared external in the calling program.  The calling 
C         sequence to MATVEC is:
C             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )
C         Where N is the number of unknowns, Y is the product A*X
C         upon return X is an input vector, NELT is the number of 
C         non-zeros in the SLAP-Column IA, JA, A storage for the matrix
C         A.  ISYM is a flag which, if non-zero, denotes that A is 
C         symmetric and only the lower or upper triangle is stored.
C MTTVEC :EXT      External.
C         Name of a routine which performs the matrix transpose vector 
C         multiply y = A'*X given A and X (where ' denotes transpose).  
C         The name of the MTTVEC routine must be declared external in
C         the calling program.  The calling sequence to MTTVEC is the
C         same as that for MATVEC, viz.:
C             CALL MTTVEC( N, X, Y, NELT, IA, JA, A, ISYM )
C         Where N is the number of unknowns, Y is the product A'*X
C         upon return X is an input vector, NELT is the number of 
C         non-zeros in the SLAP-Column IA, JA, A storage for the matrix 
C         A.  ISYM is a flag which, if non-zero, denotes that A is 
C         symmetric and only the lower or upper triangle is stored.
C MSOLVE :EXT      External.
C         Name of a routine which solves a linear system MZ = R for
C         Z given R with the preconditioning matrix M (M is supplied via
C         RWORK and IWORK arrays).  The name of the MSOLVE routine must 
C         be declared external in the calling program.  The calling 
C         sequence to MSOLVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         Where N is the number of unknowns, R is the right-hand side 
C         vector, and Z is the solution upon return.  RWORK is a real 
C         array that can be used to pass necessary preconditioning 
C         information and/or workspace to MSOLVE.  IWORK is an integer 
C         work array for the same purpose as RWORK.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual 
C         divided by the 2-norm of the right-hand side is less than TOL.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
C         residual divided by the 2-norm of M-inv times the right hand 
C         side is less than TOL, where M-inv is the inverse of the 
C         diagonal of A.
C         ITOL=11 is often useful for checking and comparing different 
C         routines.  For this case, the user must supply the "exact" 
C         solution or a very accurate approximation (one with an error 
C         much less than TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C         if ITOL=11, iteration stops when the 2-norm of the difference 
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution 
C         is less than TOL.  Note that this requires the user to set up
C         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine. 
C         The routine with this declaration should be loaded before the
C         stop test so that the correct length is used by the loader.  
C         This procedure is not standard Fortran and may not work 
C         correctly on your system (although it has worked on every
C         system the authors have tried).  If ITOL is not 11 then this
C         common block is indeed standard Fortran.
C TOL    :IN       Real.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX+1 if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Real.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Return error flag.
C           IERR = 0 => All went well.
C           IERR = 1 => Insufficient storage allocated 
C                       for WORK or IWORK.
C           IERR = 2 => Method failed to converge in 
C                       ITMAX steps.
C           IERR = 3 => Error in user input.  Check input
C                       value of N, ITOL.
C           IERR = 4 => User error tolerance set too tight.
C                       Reset to 500.0*R1MACH(3).  Iteration proceeded.
C           IERR = 5 => Preconditioning matrix, M,  is not 
C                       Positive Definite.  $(r,z) < 0.0$.
C           IERR = 6 => Matrix A is not Positive Definite.  
C                       $(p,Ap) < 0.0$.
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C R      :WORK     Real R(N).
C Z      :WORK     Real Z(N).
C P      :WORK     Real P(N).
C ATP    :WORK     Real ATP(N).
C ATZ    :WORK     Real ATZ(N).
C DZ     :WORK     Real DZ(N).
C ATDZ   :WORK     Real ATDZ(N).
C RWORK  :WORK     Real RWORK(USER DEFINABLE).
C         Real array that can be used by  MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINABLE).
C         Integer array that can be used by  MSOLVE.
C
C *Description:
C       This  routine applies the  preconditioned conjugate gradient
C       (PCG) method to a non-symmetric system of equations Ax=b. To
C       do this the normal equations are solved:
C               AA' y  = b, where  x  = A'y.
C       In PCG method the iteration count is determined by condition
C                               -1
C       number of the  matrix (M  A).   In the  situation where  the
C       normal equations are  used  to solve a  non-symmetric system
C       the condition number depends on  AA' and should therefore be
C       much worse than that of A.  This is the conventional wisdom.
C       When one has a good preconditioner for AA' this may not hold.
C       The latter is the situation when SCGN should be tried.
C       
C       If one is trying to solve  a symmetric system, SCG should be
C       used instead.
C
C       This routine does  not care  what matrix data   structure is
C       used for  A and M.  It simply   calls  the MATVEC and MSOLVE
C       routines, with  the arguments as  described above.  The user
C       could write any type of structure and the appropriate MATVEC
C       and MSOLVE routines.  It is assumed  that A is stored in the
C       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
C       stored  in  IWORK  and  RWORK)  in  some fashion.   The SLAP
C       routines SSDCGN and SSLUCN are examples of this procedure.
C       
C       Two  examples  of  matrix  data structures  are the: 1) SLAP
C       Triad  format and 2) SLAP Column format.
C       
C       =================== S L A P Triad format ===================
C
C       In  this   format only the  non-zeros are  stored.  They may
C       appear  in *ANY* order.   The user  supplies three arrays of
C       length NELT, where  NELT  is the number  of non-zeros in the
C       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero
C       the  user puts   the row  and  column index   of that matrix
C       element in the IA and JA arrays.  The  value of the non-zero
C       matrix  element is  placed in  the corresponding location of
C       the A  array.  This is  an extremely easy data  structure to
C       generate.  On  the other hand it  is  not too  efficient  on
C       vector  computers   for the  iterative  solution  of  linear
C       systems.  Hence, SLAP  changes this input  data structure to
C       the SLAP   Column  format for the  iteration (but   does not
C       change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       =================== S L A P Column format ==================
C
C       In  this format   the non-zeros are    stored counting  down
C       columns (except  for the diagonal  entry, which must  appear
C       first in each "column") and are  stored in the real array A.
C       In other words,  for  each column    in the matrix   put the
C       diagonal  entry  in A.   Then   put  in the  other  non-zero
C       elements going   down the  column (except  the  diagonal) in
C       order.  The IA array holds the row index  for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning   of   each  column.      That is,   IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1), A(JA(ICOL+1)-1)  points to the
C       end of the ICOL-th column.  Note that we always have JA(N+1)
C       = NELT+1, where N is the number of columns in the matrix and
C       NELT is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       
C *Precision:           Single Precision
C *See Also:
C         SSDCGN, SSLUCN, ISSCGN
C***REFERENCES  (NONE)
C***ROUTINES CALLED  MATVEC, MTTVEC, MSOLVE, ISSCGN, 
C                    SCOPY, SDOT, SAXPY, R1MACH
C***END PROLOGUE  SCGN
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      INTEGER IUNIT, IWORK(*)
      REAL    B(N), X(N), A(N), R(N), Z(N), P(N), ATP(N), ATZ(N)
      REAL    DZ(N), ATDZ(N), RWORK(*)
      EXTERNAL MATVEC, MTTVEC, MSOLVE
C         
C         Check user input.
C***FIRST EXECUTABLE STATEMENT  SCGN
      ITER = 0
      IERR = 0
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      TOLMIN = 500.0*R1MACH(3)
      IF( TOL.LT.TOLMIN ) THEN
         TOL = TOLMIN
         IERR = 4
      ENDIF
C         Calculate initial residual and pseudo-residual, and check
C         stopping criterion.
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      DO 10 I = 1, N
         R(I) = B(I) - R(I)
 10   CONTINUE
      CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
      CALL MTTVEC(N, Z, ATZ, NELT, IA, JA, A, ISYM)
C         
      IF( ISSCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, MSOLVE,
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, ATP, ATZ,
     $     DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 )
     $     GO TO 200
      IF( IERR.NE.0 ) RETURN
C         
C         ***** iteration loop *****
C         
      DO 100 K=1,ITMAX
         ITER = K
C         
C         Calculate coefficient BK and direction vector P.
         BKNUM = SDOT(N, Z, 1, R, 1)
         IF( BKNUM.LE.0.0 ) THEN
            IERR = 6
            RETURN
         ENDIF
         IF(ITER .EQ. 1) THEN
            CALL SCOPY(N, Z, 1, P, 1)
         ELSE
            BK = BKNUM/BKDEN
            DO 20 I = 1, N
               P(I) = Z(I) + BK*P(I)
 20         CONTINUE
         ENDIF
         BKDEN = BKNUM
C         
C         Calculate coefficient AK, new iterate X, new residual R,
C         and new pseudo-residual ATZ.
         IF(ITER .NE. 1) CALL SAXPY(N, BK, ATP, 1, ATZ, 1)
         CALL SCOPY(N, ATZ, 1, ATP, 1)
         AKDEN = SDOT(N, ATP, 1, ATP, 1)
         IF( AKDEN.LE.0.0 ) THEN
            IERR = 6
            RETURN
         ENDIF
         AK = BKNUM/AKDEN
         CALL SAXPY(N, AK, ATP, 1, X, 1)
         CALL MATVEC(N, ATP, Z, NELT, IA, JA, A, ISYM)
         CALL SAXPY(N, -AK, Z, 1, R, 1)
         CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
         CALL MTTVEC(N, Z, ATZ, NELT, IA, JA, A, ISYM)
C         
C         check stopping criterion.
         IF( ISSCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC,
     $        MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, 
     $        Z, P, ATP, ATZ, DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM, 
     $        SOLNRM) .NE. 0) GOTO 200
C         
 100  CONTINUE
C         
C         *****   end of loop  *****
C         
C         stopping criterion not satisfied.
      ITER = ITMAX + 1
C         
 200  RETURN
C------------- LAST LINE OF SCGN FOLLOWS ----------------------------
      END
*DECK SSDCGN
      SUBROUTINE SSDCGN(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
     $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
C***BEGIN PROLOGUE  SSDCGN
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSDCGN-S),
C             Non-Symmetric Linear system solve, Sparse, 
C             Iterative Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Diagonally Scaled CG Sparse Ax=b Solver for Normal Eqn's.
C            Routine to solve a general linear system Ax = b using
C            diagonal scaling with the Conjugate  Gradient  method 
C            applied to the the normal equations, viz.,  AA'y = b, 
C            where x = A'y.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C     INTEGER ITER, IERR, IUNIT, LENW, IWORK, LENIW
C     REAL B(N), X(N), A(NELT), TOL, ERR, RWORK(8*N)
C
C     CALL SSDCGN(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
C    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW)
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right-hand side vector.
C X      :INOUT    Real X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Integer A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See "Description", 
C         below.  If the SLAP Triad format is chosen it is changed 
C         internally to the SLAP Column format.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual 
C         divided by the 2-norm of the right-hand side is less than TOL.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
C         residual divided by the 2-norm of M-inv times the right hand 
C         side is less than TOL, where M-inv is the inverse of the 
C         diagonal of A.
C         ITOL=11 is often useful for checking and comparing different 
C         routines.  For this case, the user must supply the "exact" 
C         solution or a very accurate approximation (one with an error 
C         much less than TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C         if ITOL=11, iteration stops when the 2-norm of the difference 
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution 
C         is less than TOL.  Note that this requires the user to set up
C         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine. 
C         The routine with this declaration should be loaded before the
C         stop test so that the correct length is used by the loader.  
C         This procedure is not standard Fortran and may not work 
C         correctly on your system (although it has worked on every
C         system the authors have tried).  If ITOL is not 11 then this
C         common block is indeed standard Fortran.
C TOL    :IN       Real.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX+1 if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Real.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Return error flag.
C           IERR = 0 => All went well.
C           IERR = 1 => Insufficient storage allocated 
C                       for WORK or IWORK.
C           IERR = 2 => Method failed to converge in 
C                       ITMAX steps.
C           IERR = 3 => Error in user input.  Check input
C                       value of N, ITOL.
C           IERR = 4 => User error tolerance set too tight.
C                       Reset to 500.0*R1MACH(3).  Iteration proceeded.
C           IERR = 5 => Preconditioning matrix, M,  is not 
C                       Positive Definite.  $(r,z) < 0.0$.
C           IERR = 6 => Matrix A is not Positive Definite.  
C                       $(p,Ap) < 0.0$.
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C RWORK  :WORK     Real RWORK(LENW).
C         Real array used for workspace.
C LENW   :IN       Integer.
C         Length of the real workspace, RWORK.  LENW >= 8*N.
C IWORK  :WORK     Integer IWORK(LENIW).
C         Used to hold pointers into the RWORK array.
C         Upon return the following locations of IWORK hold information
C         which may be of use to the user:
C         IWORK(9)  Amount of Integer workspace actually used.
C         IWORK(10) Amount of Real    workspace actually used.
C LENIW  :IN       Integer.
C         Length of the integer workspace, IWORK.  LENIW >= 10.
C         
C *Description:
C       This  routine is simply a driver  for the  SCGN routine.  It
C       calls the   SSD2S  routine to set up the preconditioning and 
C       then calls SCGN with the appropriate   MATVEC  and    MSOLVE  
C       routines.
C
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
C       Column format.  The user can hand this routine either of the
C       of these data structures and SLAP  will figure out  which on
C       is being used and act accordingly.
C       
C       =================== S L A P Triad format ===================
C
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero   matrix  element is  placed  in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix      SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C       =================== S L A P Column format ==================
C
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C *Precision:           Single Precision
C *Side Effects:
C       The SLAP Triad format (IA, JA, A) is modified internally to be
C       the SLAP Column format.  See above.
C       
C *See Also:
C       SCGN, SSD2S, SSMV, SSMTV, SSDI
C***REFERENCES  (NONE)
C***ROUTINES CALLED  SS2Y, SCHKW, SSD2S, SCGN, SSMV, SSMTV, SSDI
C***END PROLOGUE  SSDCGN
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL
      INTEGER ITMAX, ITER, IERR, IUNIT, LENW, IWORK(LENIW), LENIW
      REAL    B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL SSMV, SSMTV, SSDI
      PARAMETER (LOCRB=1, LOCIB=11)
C
C         Modify the SLAP matrix data structure to YSMP-Column.
C***FIRST EXECUTABLE STATEMENT  SSDCGN
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL SS2Y( N, NELT, IA, JA, A, ISYM )
C
C         Set up the work arrays.
C         Compute the inverse of the diagonal of AA'.  This will be
C         used as the preconditioner.
      LOCIW = LOCIB
C
      LOCD = LOCRB
      LOCR = LOCD + N
      LOCZ = LOCR + N
      LOCP = LOCZ + N
      LOCATP = LOCP + N
      LOCATZ = LOCATP + N
      LOCDZ = LOCATZ + N
      LOCATD = LOCDZ + N
      LOCW = LOCATD + N
C
C         Check the workspace allocations.
      CALL SCHKW( 'SSDCGN', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      IWORK(4) = LOCD
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
      CALL SSD2S(N, NELT, IA, JA, A, ISYM, RWORK(1))
C         
C         Perform Conjugate Gradient algorithm on the normal equations.
      CALL SCGN( N, B, X, NELT, IA, JA, A, ISYM, SSMV, SSMTV, SSDI, 
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK(LOCR), 
     $     RWORK(LOCZ), RWORK(LOCP), RWORK(LOCATP), RWORK(LOCATZ),
     $     RWORK(LOCDZ), RWORK(LOCATD), RWORK, IWORK )
C         
      IF( ITER.GT.ITMAX ) IERR = 2
      RETURN
C------------- LAST LINE OF SSDCGN FOLLOWS ----------------------------
      END
*DECK SSLUCN
      SUBROUTINE SSLUCN(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
     $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
C***BEGIN PROLOGUE  SSLUCN
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(SSLUCN-S),
C             Non-Symmetric Linear system, Sparse, 
C             Iterative Incomplete LU Precondition
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Incomplete LU CG Sparse Ax=b Solver for Normal Equations.
C            Routine to solve  a general linear system Ax = b using the
C            incomplete  LU decomposition  with  the Conjugate Gradient
C            method  applied to the normal equations,  viz., AA'y =  b,
C            x=A'y.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C     INTEGER ITER, IERR, IUNIT, LENW, IWORK(NEL+NU+4*N+2), LENIW
C     REAL B(N), X(N), A(NELT), TOL, ERR, RWORK(NEL+NU+8*N)
C
C     CALL SSLUCN(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
C    $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right-hand side vector.
C X      :INOUT    Real X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Integer A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See "Description", 
C         below.  If the SLAP Triad format is chosen it is changed 
C         internally to the SLAP Column format.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual 
C         divided by the 2-norm of the right-hand side is less than TOL.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
C         residual divided by the 2-norm of M-inv times the right hand 
C         side is less than TOL, where M-inv is the inverse of the 
C         diagonal of A.
C         ITOL=11 is often useful for checking and comparing different 
C         routines.  For this case, the user must supply the "exact" 
C         solution or a very accurate approximation (one with an error 
C         much less than TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C         if ITOL=11, iteration stops when the 2-norm of the difference 
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution 
C         is less than TOL.  Note that this requires the user to set up
C         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine. 
C         The routine with this declaration should be loaded before the
C         stop test so that the correct length is used by the loader.  
C         This procedure is not standard Fortran and may not work 
C         correctly on your system (although it has worked on every
C         system the authors have tried).  If ITOL is not 11 then this
C         common block is indeed standard Fortran.
C TOL    :IN       Real.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX+1 if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Real.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Return error flag.
C           IERR = 0 => All went well.
C           IERR = 1 => Insufficient storage allocated 
C                       for WORK or IWORK.
C           IERR = 2 => Method failed to converge in 
C                       ITMAX steps.
C           IERR = 3 => Error in user input.  Check input
C                       value of N, ITOL.
C           IERR = 4 => User error tolerance set too tight.
C                       Reset to 500.0*R1MACH(3).  Iteration proceeded.
C           IERR = 5 => Preconditioning matrix, M,  is not 
C                       Positive Definite.  $(r,z) < 0.0$.
C           IERR = 6 => Matrix A is not Positive Definite.  
C                       $(p,Ap) < 0.0$.
C           IERR = 7 => Incomplete factorization broke down
C                       and was fudged.  Resulting preconditioning may
C                       be less than the best.
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C RWORK  :WORK     Real RWORK(LENW).
C         Real array used for workspace.  NEL is the number of non-
C         zeros in the lower triangle of the matrix (including the
C         diagonal).  NU is the number of nonzeros in the upper
C         triangle of the matrix (including the diagonal).
C LENW   :IN       Integer.
C         Length of the real workspace, RWORK.  LENW >= NEL+NU+8*N.
C IWORK  :WORK     Integer IWORK(LENIW).
C         Integer array used for workspace.  NEL is the number of non-
C         zeros in the lower triangle of the matrix (including the
C         diagonal).  NU is the number of nonzeros in the upper
C         triangle of the matrix (including the diagonal).
C         Upon return the following locations of IWORK hold information
C         which may be of use to the user:
C         IWORK(9)  Amount of Integer workspace actually used.
C         IWORK(10) Amount of Real    workspace actually used.
C LENIW  :IN       Integer.
C         Length of the integer workspace, IWORK.  LENIW >= 
C         NEL+NU+4*N+12.
C
C *Description:
C       This  routine is simply a driver  for the  SCGN  routine.    It
C       calls the SSILUS routine to set up the preconditioning and then
C       calls SCGN with the appropriate  MATVEC  and  MSOLVE  routines.
C       
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
C       Column format.  The user can hand this routine either of the
C       of these data structures and SLAP  will figure out  which on
C       is being used and act accordingly.
C       
C       =================== S L A P Triad format ===================
C
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero   matrix  element is  placed  in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C       =================== S L A P Column format ==================
C
C       This routine  requires that  the matrix A  be stored in  the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except  for the diagonal entry, which
C       must appear first in each  "column") and  are stored  in the
C       real array A.  In other words, for each column in the matrix
C       put the diagonal entry in A.  Then put in the other non-zero
C       elements going down   the  column (except  the diagonal)  in
C       order.  The IA array holds the row  index for each non-zero.
C       The JA array holds the offsets into the IA, A arrays for the
C       beginning of   each    column.    That  is,    IA(JA(ICOL)),
C       A(JA(ICOL)) points to the beginning of the ICOL-th column in
C       IA and  A.  IA(JA(ICOL+1)-1),  A(JA(ICOL+1)-1) points to the
C       end  of   the ICOL-th  column.  Note   that  we  always have
C       JA(N+1) = NELT+1, where  N  is the number of columns in  the
C       matrix and  NELT   is the number of non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C
C           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C *Precision:           Single Precision
C *Side Effects:
C       The SLAP Triad format (IA, JA, A) is modified internally to be
C       the SLAP Column format.  See above.
C
C *See Also:
C       SCGN, SDCGN, SSILUS
C***REFERENCES  (NONE)
C***ROUTINES CALLED  SS2Y, SSILUS, SCHKW, SSMV, SSMTV, SSMMTI, SCGN
C***END PROLOGUE  SSLUCN
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      INTEGER IERR, IUNIT, LENW, IWORK(LENIW), LENIW
      REAL    B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      PARAMETER (LOCRB=1, LOCIB=11)
C
      EXTERNAL SSMV, SSMTV, SSMMTI
C
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.
C***FIRST EXECUTABLE STATEMENT  SSLUCN
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL SS2Y( N, NELT, IA, JA, A, ISYM )
C
C         Count number of Non-Zero elements preconditioner ILU matrix.
C         Then set up the work arrays.
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
C         Don't count diagional.
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
CVD$ NOVECTOR
            DO 10 J = JBGN, JEND
               IF( IA(J).GT.ICOL ) THEN
                  NL = NL + 1
                  IF( ISYM.NE.0 ) NU = NU + 1
               ELSE
                  NU = NU + 1
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
C         
      LOCIL = LOCIB
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
C
      LOCL = LOCRB
      LOCDIN = LOCL + NL
      LOCU = LOCDIN + N
      LOCR = LOCU + NU
      LOCZ = LOCR + N
      LOCP = LOCZ + N
      LOCATP = LOCP + N
      LOCATZ = LOCATP + N
      LOCDZ = LOCATZ + N
      LOCATD = LOCDZ + N
      LOCW = LOCATD + N
C
C         Check the workspace allocations.
      CALL SCHKW( 'SSLUCN', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
C         Compute the Incomplete LU decomposition.
      CALL SSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCU), IWORK(LOCNR), IWORK(LOCNC) )
C         
C         Perform Conjugate Gradient algorithm on the normal equations.
      CALL SCGN(N, B, X, NELT, IA, JA, A, ISYM, SSMV, SSMTV, SSMMTI, 
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, RWORK(LOCR),
     $     RWORK(LOCZ), RWORK(LOCP), RWORK(LOCATP), RWORK(LOCATZ),
     $     RWORK(LOCDZ), RWORK(LOCATD), RWORK, IWORK ) 
C         
      IF( ITER.GT.ITMAX ) IERR = 2         
      RETURN
C------------- LAST LINE OF SSLUCN FOLLOWS ----------------------------
      END
*DECK ISSCGN
      FUNCTION ISSCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC,
     $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, 
     $     P, ATP, ATZ, DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM, SOLNRM)
C***BEGIN PROLOGUE  ISSCGN
C***REFER TO  SCGN, SSDCGN, SSLUCN
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=SINGLE PRECISION(ISSCGN-S),
C             Non-Symmetric Linear system, Sparse, 
C             Iterative Precondition, Normal Equations
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Preconditioned CG on Normal Equations Stop Test.
C            This routine calculates the stop test for the Conjugate
C            Gradient iteration   scheme applied  to     the  normal
C            equations.  It returns a nonzero  if the error estimate
C            (the type of which is determined by  ITOL) is less than
C            the user specified tolerance TOL.
C***DESCRIPTION
C *Usage:
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
C     INTEGER  IERR, IUNIT, IWORK(USER DEFINED)
C     REAL     B(N), X(N), A(N), TOL, ERR, R(N), Z(N), P(N), ATP(N)
C     REAL     ATZ(N), DZ(N), ATDZ(N)
C     REAL     RWORK(USER DEFINED), AK, BK, BNRM, SOLNRM
C     EXTERNAL MATVEC, MTTVEC, MSOLVE
C
C     IF( ISTPCGN(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MTTVEC, 
C    $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, Z, P, 
C    $     ATP, ATZ, DZ, ATDZ, RWORK, IWORK, AK, BK, BNRM, SOLNRM) 
C    $     .NE. 0 ) THEN ITERATION DONE
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Real B(N).
C         Right-hand side vector.
C X      :IN       Real X(N).
C         The current approximate solution vector.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Real A(NELT).
C         These arrays contain the matrix data structure for A.
C         It could take any form.  See "Description" in the 
C         SSCGN routine.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MATVEC :EXT      External.
C         Name of a routine which performs the matrix vector multiply
C         Y = A*X given A and X.  The name of the MATVEC routine must 
C         be declared external in the calling program.  The calling 
C         sequence to MATVEC is:
C             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )
C         Where N is the number of unknowns, Y is the product A*X
C         upon return X is an input vector, NELT is the number of 
C         non-zeros in the SLAP-Column IA, JA, A storage for the matrix 
C         A.  ISYM is a flag which, if non-zero, denotes that A is 
C         symmetric and only the lower or upper triangle is stored.
C MTTVEC :EXT      External.
C         Name of a routine which performs the matrix transpose vector 
C         multiply y = A'*X given A and X (where ' denotes transpose).  
C         The name of the MTTVEC routine must be declared external in
C         the calling program.  The calling sequence to MTTVEC is the
C         same as that for MATVEC, viz.:
C             CALL MTTVEC( N, X, Y, NELT, IA, JA, A, ISYM )
C         Where N is the number of unknowns, Y is the product A'*X
C         upon return X is an input vector, NELT is the number of 
C         non-zeros in the SLAP-Column IA, JA, A storage for the matrix 
C         A.  ISYM is a flag which, if non-zero, denotes that A is 
C         symmetric and only the lower or upper triangle is stored.
C MSOLVE :EXT      External.
C         Name of a routine which solves a linear system MZ = R for
C         Z given R with the preconditioning matrix M (M is supplied via
C         RWORK and IWORK arrays).  The name of the MSOLVE routine must 
C         be declared external in the calling program.  The calling 
C         sequence to MSOLVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         Where N is the number of unknowns, R is the right-hand side 
C         vector, and Z is the solution upon return.  RWORK is a real 
C         array that can be used to pass necessary preconditioning 
C         information and/or workspace to MSOLVE.  IWORK is an integer 
C         work array for the same purpose as RWORK.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual 
C         divided by the 2-norm of the right-hand side is less than TOL.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
C         residual divided by the 2-norm of M-inv times the right hand 
C         side is less than tol, where M-inv is the inverse of the 
C         diagonal of A.
C         ITOL=11 is often useful for checking and comparing different 
C         routines.  For this case, the user must supply the "exact" 
C         solution or a very accurate approximation (one with an error 
C         much less than TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C         if ITOL=11, iteration stops when the 2-norm of the difference 
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution 
C         is less than TOL.  Note that this requires the user to set up
C         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine. 
C         The routine with this declaration should be loaded before the
C         stop test so that the correct length is used by the loader.  
C         This procedure is not standard Fortran and may not work 
C         correctly on your system (although it has worked on every
C         system the authors have tried).  If ITOL is not 11 then this
C         common block is indeed standard Fortran.
C TOL    :IN       Real.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :IN       Integer.
C         The iteration for which to check for convergence.
C ERR    :OUT      Real.
C         Error estimate of error in the X(N) approximate solution, as 
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Error flag.  IERR is set to 3 if ITOL is not on of the 
C         acceptable values, see above. 
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C R      :IN       Real R(N).
C         The residual R = B-AX.
C Z      :WORK     Real Z(N).
C P      :IN       Real P(N).
C         The conjugate direction vector.
C ATP    :IN       Real ATP(N).
C         A-transpose times the conjugate direction vector.
C ATZ    :IN       Real ATZ(N).
C         A-transpose times the pseudo-residual.
C DZ     :IN       Real DZ(N).
C         Workspace used to hold temporary vector(s).
C ATDZ   :WORK       Real ATDZ(N).
C         Workspace.
C RWORK  :WORK     Real RWORK(USER DEFINABLE).
C         Real array that can be used by MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINABLE).
C         Integer array that can be used by MSOLVE.
C BNRM   :INOUT    Real.
C         Norm of the right hand side.  Type of norm depends on ITOL.
C         Calculated only on the first call.
C SOLNRM :INOUT    Real.
C         2-Norm of the true solution, SOLN.  Only computed and used
C         if ITOL = 11.
C
C *Function Return Values:
C       0 : Error estimate (determined by ITOL) is *NOT* less than the 
C           specified tolerance, TOL.  The iteration must continue.
C       1 : Error estimate (determined by ITOL) is less than the 
C           specified tolerance, TOL.  The iteration can be considered
C           complete.
C
C *Precision:           Single Precision
C *See Also:
C       SSCGN
C
C *Cautions:
C     This routine will attempt to write to the fortran logical output 
C     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
C     this  logical  unit  must  be  attached  to  a  file or terminal
C     before calling this routine with a non-zero  value  for   IUNIT.
C     This routine does not check for the validity of a non-zero IUNIT
C     unit number.
C***REFERENCES  (NONE)
C***ROUTINES CALLED  MATVEC, MTTVEC, MSOLVE and the BLAS
C***COMMON BLOCKS    SOLBLK
C***END PROLOGUE  ISSCGN
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      INTEGER IUNIT, IWORK(*)
      REAL    B(N), X(N), A(N), TOL, ERR, R(N), Z(N), P(N), ATP(N)
      REAL    ATZ(N), DZ(N), ATDZ(N), RWORK(*), AK, BK, BNRM, SOLNRM
      EXTERNAL MATVEC, MTTVEC, MSOLVE
      COMMON /SOLBLK/ SOLN(1)
C         
C***FIRST EXECUTABLE STATEMENT  ISSCGN
      ISSCGN = 0
C     
      IF( ITOL.EQ.1 ) THEN
C         err = ||Residual||/||RightHandSide|| (2-Norms).
         IF(ITER .EQ. 0) BNRM = SNRM2(N, B, 1)
         ERR = SNRM2(N, R, 1)/BNRM
      ELSE IF( ITOL.EQ.2 ) THEN
C                  -1              -1
C         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).
         IF(ITER .EQ. 0) THEN
            CALL MSOLVE(N, B, DZ, NELT, IA, JA, A, ISYM, RWORK, IWORK)
            CALL MTTVEC(N, DZ, ATDZ, NELT, IA, JA, A, ISYM)
            BNRM = SNRM2(N, ATDZ, 1)
         ENDIF
         ERR = SNRM2(N, ATZ, 1)/BNRM
      ELSE IF( ITOL.EQ.11 ) THEN
C         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).
         IF(ITER .EQ. 0) SOLNRM = SNRM2(N, SOLN, 1)
         DO 10 I = 1, N
            DZ(I) = X(I) - SOLN(I)
 10      CONTINUE
         ERR = SNRM2(N, DZ, 1)/SOLNRM
      ELSE
C
C         If we get here ITOL is not one of the acceptable values.
         ERR = 1.0E10
         IERR = 3
      ENDIF
C         
      IF( IUNIT.NE.0 ) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL
         ENDIF
         WRITE(IUNIT,1010) ITER, ERR, AK, BK
      ENDIF
      IF( ERR.LE.TOL ) ISSCGN = 1
C         
      RETURN
 1000 FORMAT(' PCG Applied to the Normal Equations for ',
     $     'N, ITOL = ',I5, I5,
     $     /' ITER','   Error Estimate','            Alpha',
     $     '             Beta')
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7,1X,E16.7)
C------------- LAST LINE OF ISSCGN FOLLOWS ----------------------------
      END
