#pragma once

// maybe we will later need to include Blas.hpp
#include "../libf2c/f2c.h"

namespace LaPackCPP {


//=================================================================================================
// Driver routines:

/** gbsv (general banded solver) computes the solution to a system of linear equations A * X = B 
for general banded (GB) matrices (simple driver).

Purpose:
gbsv computes the solution to a real system of linear equations A * X = B, where A is a band matrix
of order N with KL subdiagonals and KU superdiagonals, and X and B are N-by-NRHS matrices. The LU 
decomposition with partial pivoting and row interchanges is used to factor A as A = L * U, where L 
is a product of permutation and unit lower triangular matrices with KL subdiagonals, and U is upper
triangular with KL+KU superdiagonals. The factored form of A is then used to solve the system of 
equations A * X = B.

Arguments:
N:    The number of linear equations, i.e., the order of the matrix A.  N >= 0. 
KL:   The number of subdiagonals within the band of A.  KL >= 0.
KU:   The number of superdiagonals within the band of A.  KU >= 0.
NRHS: The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
AB:   array, dimension (LDAB,N). On entry, the matrix A in band storage, in rows KL+1 to 
      2*KL+KU+1; rows 1 to KL of the array need not be set. The j-th column of A is stored in the 
      j-th column of the array AB as follows: AB(KL+KU+1+i-j,j) = A(i,j) for 
      max(1,j-KU)<=i<=min(N,j+KL). On exit, details of the factorization: U is stored as an upper 
      triangular band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and the multipliers 
      used during the factorization are stored in rows KL+KU+2 to 2*KL+KU+1. See below for further
      details.
LDAB: The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
IPIV: array, dimension (N). The pivot indices that define the permutation matrix P; row i of the 
      matrix was interchanged with row IPIV(i).
B:    array, dimension (LDB,NRHS). On entry, the N-by-NRHS right hand side matrix B. On exit, if 
      INFO = 0, the N-by-NRHS solution matrix X.
LDB:  The leading dimension of the array B.  LDB >= max(1,N).
INFO  = 0: successful exit
      < 0: if INFO = -i, the i-th argument had an illegal value
      > 0: if INFO = i, U(i,i) is exactly zero.  The factorization has been completed, but the 
           factor U is exactly singular, and the solution has not been computed.

Further Details: 
The band storage scheme is illustrated by the following example, when M = N = 6, KL = 2, KU = 1:

On entry:                        On exit: 

 *    *    *    +    +    +       *    *    *   u14  u25  u36 
 *    *    +    +    +    +       *    *   u13  u24  u35  u46 
 *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56 
a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66 
a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   * 
a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    * 

Array elements marked * are not used by the routine; elements marked + need not be set on entry, 
but are required by the routine to store elements of U because of fill-in resulting from the row 
interchanges. */
template<class T>
int gbsv(long int *n, long int *kl, long int *ku, long int *nrhs, T *ab, long int *ldab, 
  long int *ipiv, T *b, long int *ldb, long int *info);
// maybe change the "long int"s back to "integer" from f2c.h - or maybe change all to "long int"
// and include f2c.h only in LaPack.cpp

// More on band-storage in lapack:
// http://www.netlib.org/lapack/lug/node124.html

// comparison with linpack and eispack:
// http://www.netlib.org/lapack/lug/node147.html

//-------------------------------------------------------------------------------------------------

/**
Purpose:
gbsvx uses the LU factorization to compute the solution to a real system of linear equations 
A * X = B, A**T * X = B, or A**H * X = B, where A is a band matrix of order N with KL subdiagonals 
and KU superdiagonals, and X and B are N-by-NRHS matrices. Error bounds on the solution and a 
condition estimate are also provided.

Description:
The following steps are performed by this subroutine:
1. If FACT = 'E', real scaling factors are computed to equilibrate the system:
      TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
      TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
      TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
   Whether or not the system will be equilibrated depends on the scaling of the matrix A, but if 
   equilibration is used, A is overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
   or diag(C)*B (if TRANS = 'T' or 'C').
2. If FACT = 'N' or 'E', the LU decomposition is used to factor the matrix A (after equilibration 
   if FACT = 'E') as
      A = L * U,
   where L is a product of permutation and unit lower triangular matrices with KL subdiagonals, and
   U is upper triangular with KL+KU superdiagonals.
3. If some U(i,i)=0, so that U is exactly singular, then the routine returns with INFO = i. 
   Otherwise, the factored form of A is used to estimate the condition number of the matrix A.  If 
   the reciprocal of the condition number is less than machine precision, INFO = N+1 is returned as
   a warning, but the routine still goes on to solve for X and compute error bounds as described below.
4. The system of equations is solved for X using the factored form of A.
5. Iterative refinement is applied to improve the computed solution matrix and calculate error 
   bounds and backward error estimates for it.
6. If equilibration was used, the matrix X is premultiplied by diag(C) (if TRANS = 'N') or diag(R) 
   (if TRANS = 'T' or 'C') so that it solves the original system before equilibration.

Arguments:
FACT:  Specifies whether or not the factored form of the matrix A is supplied on entry, and if not, 
       whether the matrix A should be equilibrated before it is factored.
       = 'F':  On entry, AFB and IPIV contain the factored form of A.  If EQUED is not 'N', the 
               matrix A has been equilibrated with scaling factors given by R and C. AB, AFB, and 
               IPIV are not modified.
       = 'N':  The matrix A will be copied to AFB and factored.
       = 'E':  The matrix A will be equilibrated if necessary, then copied to AFB and factored.
TRANS: Specifies the form of the system of equations.
       = 'N':  A * X = B     (No transpose)
       = 'T':  A**T * X = B  (Transpose)
       = 'C':  A**H * X = B  (Transpose) ...ToDo: enable for complex use -> hermitian transpose
N:     The number of linear equations, i.e., the order of the matrix A.  N >= 0.
KL:    The number of subdiagonals within the band of A.  KL >= 0.
KU:    The number of superdiagonals within the band of A.  KU >= 0.
NRHS:  The number of right hand sides, i.e., the number of columns of the matrices B and X.  NRHS >= 0.
AB:    array, dimension (LDAB,N). On entry, the matrix A in band storage, in rows 1 to KL+KU+1. The 
       j-th column of A is stored in the j-th column of the array AB as follows:
       AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl). If FACT = 'F' and EQUED is not 'N', 
       then A must have been equilibrated by the scaling factors in R and/or C.  AB is not modified
       if FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit. On exit, if EQUED .ne. 'N', 
       A is scaled as follows:
       EQUED = 'R':  A := diag(R) * A
       EQUED = 'C':  A := A * diag(C)
       EQUED = 'B':  A := diag(R) * A * diag(C).
LDAB:  The leading dimension of the array AB.  LDAB >= KL+KU+1.
AFB:   array, dimension (LDAFB,N). If FACT = 'F', then AFB is an input argument and on entry 
       contains details of the LU factorization of the band matrix A, as computed by DGBTRF. U is 
       stored as an upper triangular band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, 
       and the multipliers used during the factorization are stored in rows KL+KU+2 to 2*KL+KU+1. 
       If EQUED .ne. 'N', then AFB is the factored form of the equilibrated matrix A. If FACT = 'N', 
       then AFB is an output argument and on exit returns details of the LU factorization of A. If 
       FACT = 'E', then AFB is an output argument and on exit returns details of the LU 
       factorization of the equilibrated matrix A (see the description of AB for the form of the
       equilibrated matrix).
LDAFB: The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.
IPIV:  array, dimension (N). If FACT = 'F', then IPIV is an input argument and on entry contains 
       the pivot indices from the factorization A = L*U as computed by DGBTRF; row i of the matrix 
       was interchanged with row IPIV(i). If FACT = 'N', then IPIV is an output argument and on 
       exit contains the pivot indices from the factorization A = L*U of the original matrix A. If 
       FACT = 'E', then IPIV is an output argument and on exit contains the pivot indices from the 
       factorization A = L*U of the equilibrated matrix A.
EQUED: Specifies the form of equilibration that was done.
       = 'N':  No equilibration (always true if FACT = 'N').
       = 'R':  Row equilibration, i.e., A has been premultiplied by diag(R).
       = 'C':  Column equilibration, i.e., A has been postmultiplied by diag(C).
       = 'B':  Both row and column equilibration, i.e., A has been replaced by diag(R) * A * diag(C).
EQUED: is an input argument if FACT = 'F'; otherwise, it is an output argument.
R:     array, dimension (N). The row scale factors for A.  If EQUED = 'R' or 'B', A is multiplied 
       on the left by diag(R); if EQUED = 'N' or 'C', R is not accessed.  R is an input argument if
       FACT = 'F'; otherwise, R is an output argument.  If FACT = 'F' and EQUED = 'R' or 'B', each 
       element of R must be positive.
C:     array, dimension (N). The column scale factors for A.  If EQUED = 'C' or 'B', A is 
       multiplied on the right by diag(C); if EQUED = 'N' or 'R', C is not accessed.  C is an input
       argument if FACT = 'F'; otherwise, C is an output argument.  If FACT = 'F' and 
       EQUED = 'C' or 'B', each element of C must be positive.
B:     array, dimension (LDB,NRHS). On entry, the right hand side matrix B. On exit, if 
       EQUED = 'N', B is not modified; if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
       diag(R)*B; if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is overwritten by diag(C)*B.
LDB:   The leading dimension of the array B.  LDB >= max(1,N).
X:     array, dimension (LDX,NRHS). If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X to 
       the original system of equations.  Note that A and B are modified on exit if EQUED .ne. 'N', 
       and the solution to the equilibrated system is inv(diag(C))*X if TRANS = 'N' and 
       EQUED = 'C' or 'B', or inv(diag(R))*X if TRANS = 'T' or 'C' and EQUED = 'R' or 'B'.
LDX:   The leading dimension of the array X.  LDX >= max(1,N).
RCOND: The estimate of the reciprocal condition number of the matrix A after equilibration
        (if done). If RCOND is less than the machine precision (in particular, if RCOND = 0), the
        matrix is singular to working precision.  This condition is indicated by a return code of 
        INFO > 0.
FERR:   array, dimension (NRHS). The estimated forward error bound for each solution vector
        X(j) (the j-th column of the solution matrix X). If XTRUE is the true solution 
        corresponding to X(j), FERR(j) is an estimated upper bound for the magnitude of the largest
        element in (X(j) - XTRUE) divided by the magnitude of the largest element in X(j). The 
        estimate is as reliable as the estimate for RCOND, and is almost always a slight 
        overestimate of the true error.
BERR:   array, dimension (NRHS). The componentwise relative backward error of each solution vector
        X(j) (i.e., the smallest relative change in any element of A or B that makes X(j) an exact 
        solution).
WORK:   array, dimension (3*N). On exit, WORK(1) contains the reciprocal pivot growth factor 
        norm(A)/norm(U). The "max absolute element" norm is used. If WORK(1) is much less than 1, 
        then the stability of the LU factorization of the (equilibrated) matrix A could be poor. 
        This also means that the solution X, condition estimator RCOND, and forward error bound 
        FERR could be unreliable. If factorization fails with 0<INFO<=N, then WORK(1) contains the 
        reciprocal pivot growth factor for the leading INFO columns of A.
IWORK:  array, dimension (N)
INFO:   = 0: successful exit
        < 0: if INFO = -i, the i-th argument had an illegal value
        > 0: if INFO = i, and i is
             <= N: U(i,i) is exactly zero. The factorization has been completed, but the factor U 
                   is exactly singular, so the solution and error bounds could not be computed. 
                   RCOND = 0 is returned.
             = N+1: U is nonsingular, but RCOND is less than machine precision, meaning that the 
                    matrix is singular to working precision.  Nevertheless, the solution and error 
                    bounds are computed because there are a number of situations where the computed
                    solution can be more accurate than the value of RCOND would suggest. */
template<class T>
int gbsvx(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, 
  doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, integer *ipiv, char *equed, 
  doublereal *r__, doublereal *c__, doublereal *b, integer *ldb, doublereal *x, integer *ldx, 
  doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, 
  integer *info, ftnlen fact_len, ftnlen trans_len, ftnlen equed_len);


//=================================================================================================
// Computational routines:



//-------------------------------------------------------------------------------------------------

/**
Purpose: 
GBTF2 computes an LU factorization of a real m-by-n band matrix A using partial pivoting 
with row interchanges. This is the unblocked version of the algorithm, calling Level 2 BLAS.

Arguments:
M:    The number of rows of the matrix A.  M >= 0.
N:    The number of columns of the matrix A.  N >= 0.
KL:   The number of subdiagonals within the band of A.  KL >= 0.
KU:   The number of superdiagonals within the band of A.  KU >= 0.
AB:   array, dimension (LDAB,N). On entry, the matrix A in band storage, in rows KL+1 to 
      2*KL+KU+1; rows 1 to KL of the array need not be set. The j-th column of A is stored in the 
      j-th column of the array AB as follows:
       AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
       On exit, details of the factorization: U is stored as an upper triangular band matrix with 
       KL+KU superdiagonals in rows 1 to KL+KU+1, and the multipliers used during the factorization 
       are stored in rows KL+KU+2 to 2*KL+KU+1. See below for further details.
LDAB:  The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
IPIV:  array, dimension (min(M,N)). The pivot indices; for 1 <= i <= min(M,N), row i of the matrix 
       was interchanged with row IPIV(i).
INFO:  = 0: successful exit
       < 0: if INFO = -i, the i-th argument had an illegal value
       > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
            has been completed, but the factor U is exactly singular, and division by zero will 
            occur if it is used to solve a system of equations.

Further Details:
The band storage scheme is illustrated by the following example, when M = N = 6, KL = 2, KU = 1:

 On entry:                         On exit:

  *    *    *    +    +    +       *    *    *   u14  u25  u36
  *    *    +    +    +    +       *    *   u13  u24  u35  u46
  *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56
 a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66
 a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   *
 a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    *

Array elements marked * are not used by the routine; elements marked + need not be set on entry, 
but are required by the routine to store elements of U, because of fill-in resulting from the 
row */
template<class T>
int gbtf2(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, integer *ipiv, 
  integer *info);

//-------------------------------------------------------------------------------------------------

/** gbtrf computes an LU factorization of a real m-by-n band matrix A.

Purpose:
gbtrf computes an LU factorization of a real m-by-n band matrix A using partial pivoting with row 
interchanges. This is the blocked version of the algorithm, calling Level 3 BLAS.

Arguments:
M:    The number of rows of the matrix A.  M >= 0.
N:    The number of columns of the matrix A.  N >= 0.
KL:   The number of subdiagonals within the band of A.  KL >= 0. 
KU:   The number of superdiagonals within the band of A.  KU >= 0.
AB:   array, dimension (LDAB,N). On entry, the matrix A in band storage, in rows KL+1 to 2*KL+KU+1; 
      rows 1 to KL of the array need not be set. The j-th column of A is stored in the j-th column 
      of the array AB as follows: AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl) 
      On exit, details of the factorization: U is stored as an upper triangular band matrix with 
      KL+KU superdiagonals in rows 1 to KL+KU+1, and the multipliers used during the factorization 
      are stored in rows KL+KU+2 to 2*KL+KU+1. See below for further details. 
LDAB: The leading dimension of the array AB. LDAB >= 2*KL+KU+1. 
IPIV: array, dimension (min(M,N)). The pivot indices; for 1 <= i <= min(M,N), row i of the matrix 
      was interchanged with row IPIV(i). 
INFO: = 0: successful exit 
      < 0: if INFO = -i, the i-th argument had an illegal value 
      > 0: if INFO = +i, U(i,i) is exactly zero. The factorization has been completed, but the 
           factor U is exactly singular, and division by zero will occur if it is used to solve a 
           system of equations.

Further Details:
The band storage scheme is illustrated by the following example, when M = N = 6, KL = 2, KU = 1: 

On entry:                        On exit: 
 *    *    *    +    +    +       *    *    *   u14  u25  u36 
 *    *    +    +    +    +       *    *   u13  u24  u35  u46 
 *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56 
a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66 
a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   * 
a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    * 

Array elements marked * are not used by the routine; elements marked + need not be set on entry, 
but are required by the routine to store elements of U because of fill-in resulting from the row 
interchanges. */
template<class T>
int gbtrf(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, 
  integer *ipiv, integer *info);
// TRF: triangular factorization?

//-------------------------------------------------------------------------------------------------

/**
Purpose:
gbtrs solves a system of linear equations A * X = B  or  A**T * X = B with a general band matrix A 
using the LU factorization computed by gbtrf.

Arguments: 
TRANS: Specifies the form of the system of equations.
       = 'N':  A * X = B  (No transpose)
       = 'T':  A**T* X = B  (Transpose) 
       = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
N:     The order of the matrix A.  N >= 0. 
KL:    The number of subdiagonals within the band of A.  KL >= 0.
KU:    The number of superdiagonals within the band of A.  KU >= 0.
NRHS:  The number of right hand sides, i.e., the number of columns of the matrix B.  NRHS >= 0.
AB:    array, dimension (LDAB,N). Details of the LU factorization of the band matrix A, as computed 
       by gbtrf. U is stored as an upper triangular band matrix with KL+KU superdiagonals in rows 1 
       to KL+KU+1, and the multipliers used during the factorization are stored in rows KL+KU+2 to 
       2*KL+KU+1.
LDAB:  The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
IPIV:  array, dimension (N). The pivot indices; for 1 <= i <= N, row i of the matrix was 
       interchanged with row IPIV(i).
B:     array, dimension (LDB,NRHS). On entry, the right hand side matrix B. On exit, the solution 
       matrix X. 
LDB:   The leading dimension of the array B.  LDB >= max(1,N).
INFO:  = 0:  successful exit
       < 0: if INFO = -i, the i-th argument had an illegal value */
template<class T>
int gbtrs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab, integer *ldab, 
  integer *ipiv, T *b, integer *ldb, integer *info, ftnlen trans_len);
// TRS: triangular solution or back(S)ubstitution? -> look up in manual

//=================================================================================================
// Auxiliary routines:




//-------------------------------------------------------------------------------------------------

/** IEEECK is called from the ILAENV to verify that Infinity and possibly NaN arithmetic is safe
(i.e. will not trap).

Arguments:
ISPEC: Specifies whether to test just for inifinity arithmetic or whether to test for infinity 
       and NaN arithmetic.
       = 0: Verify infinity arithmetic only.
       = 1: Verify infinity and NaN arithmetic.
ZERO:  Must contain the value 0.0. This is passed to prevent the compiler from optimizing away this
       code.
ONE:   Must contain the value 1.0. This is passed to prevent the compiler from optimizing away this 
       code.
RETURN VALUE: = 0:  Arithmetic failed to produce the correct answers
              = 1:  Arithmetic produced the correct answers  */
integer ieeeck(integer *ispec, f2c_real *zero, f2c_real *one);

//-------------------------------------------------------------------------------------------------

/**

Purpose:
ILAENV is called from the LAPACK routines to choose problem-dependent parameters for the local 
environment.  See ISPEC for a description of the parameters.

ILAENV returns an INTEGER
if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC
if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value.

This version provides a set of parameters which should give good, but not optimal, performance on
many of the currently available computers.  Users are encouraged to modify this subroutine to set
the tuning parameters for their particular machine using the option and problem size information in
the arguments.

This routine will not function correctly if it is converted to all lower case. Converting it to all
upper case is allowed.

Arguments:
ISPEC: Specifies the parameter to be returned as the value of ILAENV.
       = 1: the optimal blocksize; if this value is 1, an unblocked
            algorithm will give the best performance.
       = 2: the minimum block size for which the block routine
            should be used; if the usable block size is less than
            this value, an unblocked routine should be used.
       = 3: the crossover point (in a block routine, for N less
            than this value, an unblocked routine should be used)
       = 4: the number of shifts, used in the nonsymmetric
            eigenvalue routines (DEPRECATED)
       = 5: the minimum column dimension for blocking to be used;
            rectangular blocks must have dimension at least k by m,
            where k is given by ILAENV(2,...) and m by ILAENV(5,...)
       = 6: the crossover point for the SVD (when reducing an m by n
            matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
            this value, a QR factorization is used first to reduce
            the matrix to a triangular form.)
       = 7: the number of processors
       = 8: the crossover point for the multishift QR method
            for nonsymmetric eigenvalue problems (DEPRECATED)
       = 9: maximum size of the subproblems at the bottom of the
            computation tree in the divide-and-conquer algorithm
            (used by xGELSD and xGESDD)
       =10: ieee NaN arithmetic can be trusted not to trap
       =11: infinity arithmetic can be trusted not to trap
        12 <= ISPEC <= 16:
            xHSEQR or related subroutines,
            see IPARMQ for detailed explanation
NAME: The name of the calling subroutine, in either upper case or lower case.
OPTS: The character options to the subroutine NAME, concatenated into a single character string. 
      For example, UPLO = 'U', TRANS = 'T', and DIAG = 'N' for a triangular routine would be 
      specified as OPTS = 'UTN'.
N1:
N2:
N3:
N4:   Problem dimensions for the subroutine NAME; these may not all be required.

Further Details:
The following conventions have been used when calling ILAENV from the
LAPACK routines:
1)  OPTS is a concatenation of all of the character options to
    subroutine NAME, in the same order that they appear in the
    argument list for NAME, even if they are not used in determining
    the value of the parameter specified by ISPEC.
2)  The problem dimensions N1, N2, N3, N4 are specified in the order
    that they appear in the argument list for NAME.  N1 is used
    first, N2 second, and so on, and unused problem dimensions are
    passed a value of -1.
3)  The parameter value returned by ILAENV is checked for validity in
    the calling subroutine.  For example, ILAENV is used to retrieve
    the optimal blocksize for STRTRI as follows:
    NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
    IF( NB.LE.1 ) NB = MAX( 1, N )  */
integer ilaenv(integer *ispec, char *name__, char *opts, integer *n1, 
  integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen opts_len);


//-------------------------------------------------------------------------------------------------

/**
Purpose:
This program sets problem and machine dependent parameters useful for xHSEQR and related 
subroutines for eigenvalue problems. It is called whenever IPARMQ is called 
with 12 <= ISPEC <= 16

Arguments:
ISPEC: specifies which tunable parameter IPARMQ should return.
       ISPEC=12: (INMIN)  Matrices of order nmin or less are sent directly to xLAHQR, the 
                  implicit double shift QR algorithm.  NMIN must be at least 11.
       ISPEC=13: (INWIN)  Size of the deflation window. This is best set greater than or equal to
                 the number of simultaneous shifts NS. Larger matrices benefit from larger 
                 deflation windows.
       ISPEC=14: (INIBL) Determines when to stop nibbling and invest in an (expensive) multi-shift 
                 QR sweep. If the aggressive early deflation subroutine finds LD converged 
                 eigenvalues from an order NW deflation window and LD.GT.(NW*NIBBLE)/100, then the 
                 next QR sweep is skipped and early deflation is applied immediately to the 
                 remaining active diagonal block.  Setting IPARMQ(ISPEC=14) = 0 causes TTQRE to 
                 skip a multi-shift QR sweep whenever early deflation finds a converged eigenvalue.
                 Setting IPARMQ(ISPEC=14) greater than or equal to 100 prevents TTQRE from skipping
                 a multi-shift QR sweep.
       ISPEC=15: (NSHFTS) The number of simultaneous shifts in a multi-shift QR iteration.
       ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the following meanings.
                 0: During the multi-shift QR/QZ sweep, blocked eigenvalue reordering, blocked
                    Hessenberg-triangular reduction, reflections and/or rotations are not 
                    accumulated when updating the far-from-diagonal matrix entries.
                 1: During the multi-shift QR/QZ sweep, blocked eigenvalue reordering, blocked
                    Hessenberg-triangular reduction, reflections and/or rotations are accumulated,
                    and matrix-matrix multiplication is used to update the far-from-diagonal 
                    matrix entries.
                 2: During the multi-shift QR/QZ sweep, blocked eigenvalue reordering, blocked
                    Hessenberg-triangular reduction, reflections and/or rotations are accumulated, 
                    and 2-by-2 block structure is exploited during matrix-matrix multiplies.
                    (If xTRMM is slower than xGEMM, then IPARMQ(ISPEC=16)=1 may be more efficient 
                    than IPARMQ(ISPEC=16)=2 despite the greater level of arithmetic work implied 
                    by the latter choice.)
NAME:  Name of the calling subroutine
OPTS:  This is a concatenation of the string arguments to TTQRE.
N:     the order of the Hessenberg matrix H.
ILO:
IHI:   It is assumed that H is already upper triangular in rows and columns 1:ILO-1 and IHI+1:N.
LWORK: The amount of workspace available.

Further Details:
Little is known about how best to choose these parameters. It is possible to use different values 
of the parameters for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.

It is probably best to choose different parameters for different matrices and different parameters
at different times during the iteration, but this has not been implemented --- yet.

The best choices of most of the parameters depend in an ill-understood way on the relative 
execution rate of xLAQR3 and xLAQR5 and on the nature of each particular eigenvalue problem. 
Experiment may be the only practical way to determine which choices are most effective.

Following is a list of default values supplied by IPARMQ. These defaults may be adjusted in order 
to attain better performance in any particular computational environment.

IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point. Default: 75. (Must be at least 11.)

IPARMQ(ISPEC=13) Recommended deflation window size. This depends on ILO, IHI and NS, the number of 
simultaneous shifts returned by IPARMQ(ISPEC=15). The default for (IHI-ILO+1).LE.500 is NS. The 
default for (IHI-ILO+1).GT.500 is 3*NS/2.

IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.

IPARMQ(ISPEC=15) Number of simultaneous shifts, NS. a multi-shift QR iteration.

If IHI-ILO+1 is ...

greater than      ...but less    ... the
or equal to ...      than        default is

        0               30       NS =   2+
       30               60       NS =   4+
       60              150       NS =  10
      150              590       NS =  **
      590             3000       NS =  64
     3000             6000       NS = 128
     6000             infinity   NS = 256

(+)  By default matrices of this order are passed to the implicit double shift routine xLAHQR. 
See IPARMQ(ISPEC=12) above. These values of NS are used only in case of a rare xLAHQR failure.

(**) The asterisks (**) indicate an ad-hoc function increasing from 10 to 64.

IPARMQ(ISPEC=16) Select structured matrix multiply. (See ISPEC=16 above for details.) Default: 3.
*/
integer iparmq(integer *ispec, char *name__, char *opts, integer *n, integer *ilo, integer *ihi, 
  integer *lwork, ftnlen name_len, ftnlen opts_len);

//-------------------------------------------------------------------------------------------------

/**
Purpose:
langb returns the value of the one norm, or the Frobenius norm, or the infinity norm, or the
element of largest absolute value of an n by n band matrix A, with kl sub-diagonals and ku
super-diagonals.
 langb = ( max(abs(A(i,j))), NORM = 'M' or 'm'
         (
         ( norm1(A),         NORM = '1', 'O' or 'o'
         (
         ( normI(A),         NORM = 'I' or 'i'
         (
         ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
where  norm1  denotes the  one norm of a matrix (maximum column sum), normI denotes the infinity
norm  of a matrix  (maximum row sum) and normF  denotes the  Frobenius norm of a matrix (square
root of sum of squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.

Arguments:
NORM:  Specifies the value to be returned in DLANGB as described above.
N:     The order of the matrix A.  N >= 0.  When N = 0, DLANGB is set to zero.
KL:    The number of sub-diagonals of the matrix A.  KL >= 0.
KU:    The number of super-diagonals of the matrix A.  KU >= 0.
AB:    array, dimension (LDAB,N). The band matrix A, stored in rows 1 to KL+KU+1. The j-th
       column of A is stored in the j-th column of the array AB as follows:
       AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).
LDAB:  The leading dimension of the array AB.  LDAB >= KL+KU+1.
WORK:  array, dimension (MAX(1,LWORK)), where LWORK >= N when NORM = 'I'; otherwise, WORK is not
       referenced. */
template<class T>
T langb(char *norm, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *work, 
  ftnlen norm_len);

//-------------------------------------------------------------------------------------------------

/** Copies all or part of a two-dimensional matrix A to another matrix B.

Arguments:
UPLO: Specifies the part of the matrix A to be copied to B.
      = 'U':      Upper triangular part
      = 'L':      Lower triangular part
      Otherwise:  All of the matrix A
M:    The number of rows of the matrix A.  M >= 0.
N:    The number of columns of the matrix A.  N >= 0.
A:    array, dimension (LDA,N). The m by n matrix A.  If UPLO = 'U', only the upper triangle or 
      trapezoid is accessed; if UPLO = 'L', only the lower triangle or trapezoid is accessed.
LDA:  The leading dimension of the array A.  LDA >= max(1,M).
B:    array, dimension (LDB,N). On exit, B = A in the locations specified by UPLO.
LDB:  The leading dimension of the array B.  LDB >= max(1,M). */
template<class T>
int lacpy(char *uplo, integer *m, integer *n, T *a, integer *lda, T *b, integer *ldb, 
  ftnlen uplo_len);

//-------------------------------------------------------------------------------------------------

/**
Purpose:
DLAQGB equilibrates a general M by N band matrix A with KL subdiagonals and KU superdiagonals 
using the row and scaling factors in the vectors R and C.

Arguments:
M:      The number of rows of the matrix A.  M >= 0.
N:      The number of columns of the matrix A.  N >= 0.
KL:     The number of subdiagonals within the band of A.  KL >= 0.
KU:     The number of superdiagonals within the band of A.  KU >= 0.
AB      array, dimension (LDAB,N). On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
        The j-th column of A is stored in the j-th column of the array AB as follows:
        AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl). On exit, the equilibrated matrix, 
        in the same storage format as A.  See EQUED for the form of the equilibrated matrix.
LDAB:   The leading dimension of the array AB.  LDA >= KL+KU+1.
R:      array, dimension (M). The row scale factors for A.
C:      array, dimension (N). The column scale factors for A.
ROWCND: Ratio of the smallest R(i) to the largest R(i).
COLCND: Ratio of the smallest C(i) to the largest C(i).
AMAX:   Absolute value of largest matrix entry.
EQUED:  Specifies the form of equilibration that was done.
        = 'N': No equilibration
        = 'R': Row equilibration, i.e., A has been premultiplied by diag(R).
        = 'C': Column equilibration, i.e., A has been postmultiplied by diag(C).
        = 'B': Both row and column equilibration, i.e., A has been replaced by 
               diag(R) * A * diag(C).

Internal Parameters:
THRESH is a threshold value used to decide if row or column scaling should be done based on the 
ratio of the row or column scaling factors.  If ROWCND < THRESH, row scaling is done, and if
COLCND < THRESH, column scaling is done. LARGE and SMALL are threshold values used to decide if 
row scaling should be done based on the absolute size of the largest matrix element. 
If AMAX > LARGE or AMAX < SMALL, row scaling is done.  */
template<class T>
int laqgb(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *r__, T *c__, 
  T *rowcnd, T *colcnd, T *amax, char *equed, ftnlen equed_len);

//-------------------------------------------------------------------------------------------------

/* laswp performs a series of row interchanges on the matrix A. One row interchange is 
initiated for each of rows K1 through K2 of A.

Arguments:
N:    The number of columns of the matrix A.
A:    array, dimension (LDA,N). On entry, the matrix of column dimension N to which the row 
      interchanges will be applied. On exit, the permuted matrix.
LDA:  The leading dimension of the array A.
K1:   The first element of IPIV for which a row interchange will be done.
K2:   (K2-K1+1) is the number of elements of IPIV for which a row interchange will be done.
IPIV: array, dimension (K1+(K2-K1)*abs(INCX)). The vector of pivot indices. Only the elements in 
      positions K1 through K1+(K2-K1)*abs(INCX) of IPIV are accessed. IPIV(K1+(K-K1)*abs(INCX)) = L 
      implies rows K and L are to be interchanged.
INCX: The increment between successive values of IPIV. If INCX is negative, the pivots are applied 
      in reverse order.  */
template<class T>
int laswp(integer *n, T *a, integer *lda, integer *k1, integer *k2, integer *ipiv, integer *incx);

// maybe the LaPack source files should be split into LaPackGB, etc. or: LaPackDrivers, 
// LapackComputation, LaPackAuxiliary


}