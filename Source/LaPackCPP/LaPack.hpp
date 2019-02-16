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

//-------------------------------------------------------------------------------------------------

/**
Purpose:
GBSVXX uses the LU factorization to compute the solution to a double precision system of linear
equations  A * X = B,  where A is an N-by-N matrix and X and B are N-by-NRHS matrices. If 
requested, both normwise and maximum componentwise error bounds are returned. DGBSVXX will return
a solution with a tiny guaranteed error (O(eps) where eps is the working machine precision) unless
the matrix is very ill-conditioned, in which case a warning is returned. Relevant condition 
numbers also are calculated and returned.
DGBSVXX accepts user-provided factorizations and equilibration factors; see the definitions of the
FACT and EQUED options. Solving with refinement and using a factorization from a previous DGBSVXX 
call will also produce a solution with either O(eps) errors or warnings, but we cannot make that 
claim for general user-provided factorizations and equilibration factors if they differ from what 
DGBSVXX would itself produce.

Description:
The following steps are performed:
1. If FACT = 'E', double precision scaling factors are computed to equilibrate
   the system:
     TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
     TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
     TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
   Whether or not the system will be equilibrated depends on the scaling of the matrix A, but if 
   equilibration is used, A is overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
   or diag(C)*B (if TRANS = 'T' or 'C').
2. If FACT = 'N' or 'E', the LU decomposition is used to factor the matrix A (after equilibration 
   if FACT = 'E') as A = P * L * U, where P is a permutation matrix, L is a unit lower triangular
   matrix, and U is upper triangular.
3. If some U(i,i)=0, so that U is exactly singular, then the routine returns with INFO = i. 
   Otherwise, the factored form of A is used to estimate the condition number of the matrix A (see
   argument RCOND). If the reciprocal of the condition number is less than machine precision, the
   routine still goes on to solve for X and compute error bounds as described below.
4. The system of equations is solved for X using the factored form of A.
5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero), the routine will use iterative 
   refinement to try to get a small error and error bounds.  Refinement calculates the residual to 
   at least twice the working precision.
6. If equilibration was used, the matrix X is premultiplied by diag(C) (if TRANS = 'N') or diag(R) 
   (if TRANS = 'T' or 'C') so that it solves the original system before equilibration.


Arguments:
Some optional parameters are bundled in the PARAMS array. These settings determine how refinement 
is performed, but often the defaults are acceptable. If the defaults are acceptable, users can pass
NPARAMS = 0 which prevents the source code from accessing the PARAMS argument.

FACT:  Specifies whether or not the factored form of the matrix A is supplied on entry, and if not, 
  whether the matrix A should be equilibrated before it is factored.
  = 'F':  On entry, AF and IPIV contain the factored form of A. If EQUED is not 'N', the matrix A 
     has been equilibrated with scaling factors given by R and C. A, AF, and IPIV are not modified.
  = 'N':  The matrix A will be copied to AF and factored.
  = 'E':  The matrix A will be equilibrated if necessary, then copied to AF and factored.

TRANS:  Specifies the form of the system of equations:
  = 'N':  A * X = B     (No transpose)
  = 'T':  A**T * X = B  (Transpose)
  = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)

N:    The number of linear equations, i.e., the order of the matrix A.  N >= 0.

KL:   The number of subdiagonals within the band of A.  KL >= 0.

KU:   The number of superdiagonals within the band of A.  KU >= 0.

NRHS: The number of right hand sides, i.e., the number of columns of the matrices B and X.  NRHS >= 0.

AB:   array, dimension (LDAB,N). On entry, the matrix A in band storage, in rows 1 to KL+KU+1. The 
  j-th column of A is stored in the j-th column of the array AB as follows:
    AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)
   If FACT = 'F' and EQUED is not 'N', then AB must have been equilibrated by the scaling factors 
   in R and/or C.  AB is not modified if FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit.
   On exit, if EQUED .ne. 'N', A is scaled as follows:
     EQUED = 'R':  A := diag(R) * A
     EQUED = 'C':  A := A * diag(C)
     EQUED = 'B':  A := diag(R) * A * diag(C).

LDAB: The leading dimension of the array AB.  LDAB >= KL+KU+1.

AFB:  array, dimension (LDAFB,N). If FACT = 'F', then AFB is an input argument and on entry 
  contains details of the LU factorization of the band matrix A, as computed by DGBTRF. U is stored as
  an upper triangular band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and the 
  multipliers used during the factorization are stored in rows KL+KU+2 to 2*KL+KU+1.  If 
  EQUED .ne. 'N', then AFB is the factored form of the equilibrated matrix A.
  If FACT = 'N', then AF is an output argument and on exit returns the factors L and U from the 
  factorization A = P*L*U of the original matrix A. If FACT = 'E', then AF is an output argument 
  and on exit returns the factors L and U from the factorization A = P*L*U of the equilibrated 
  matrix A (see the description of A for the form of the equilibrated matrix).

LDAFB: The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.

IPIV:  array, dimension (N). If FACT = 'F', then IPIV is an input argument and on entry contains 
  the pivot indices from the factorization A = P*L*U as computed by DGETRF; row i of the  matrix was 
  interchanged  with row IPIV(i). If FACT = 'N', then IPIV is an output argument and on exit 
  contains the pivot indices from the factorization A = P*L*U of the original matrix A. If 
  FACT = 'E', then IPIV is an output argument and on exit contains the pivot indices from the 
  factorization A = P*L*U of the equilibrated matrix A.

EQUED: Specifies the form of equilibration that was done.
  = 'N':  No equilibration (always true if FACT = 'N').
  = 'R':  Row equilibration, i.e., A has been premultiplied by diag(R).
  = 'C':  Column equilibration, i.e., A has been postmultiplied by diag(C).
  = 'B':  Both row and column equilibration, i.e., A has been replaced by diag(R) * A * diag(C).
  EQUED is an input argument if FACT = 'F'; otherwise, it is an output argument.

R:  array, dimension (N) The row scale factors for A.  If EQUED = 'R' or 'B', A is multiplied on 
  the left by diag(R); if EQUED = 'N' or 'C', R is not accessed.  R is an input argument if 
  FACT = 'F'; otherwise, R is an output argument.  If FACT = 'F' and EQUED = 'R' or 'B', each 
  element of R must be positive. If R is output, each element of R is a power of the radix. If R is
  input, each element of R should be a power of the radix to ensure a reliable solution and error 
  estimates. Scaling by powers of the radix does not cause rounding errors unless the result 
  underflows or overflows. Rounding errors during scaling lead to refining with a matrix that is 
  not equivalent to the input matrix, producing error estimates that may not be reliable.

C:  array, dimension (N). The column scale factors for A.  If EQUED = 'C' or 'B', A is multiplied 
  on the right by diag(C); if EQUED = 'N' or 'R', C is not accessed.  C is an input argument if 
  FACT = 'F'; otherwise, C is an output argument.  If FACT = 'F' and EQUED = 'C' or 'B', each 
  element of C must be positive. If C is output, each element of C is a power of the radix. If C 
  is input, each element of C should be a power of the radix to ensure a reliable solution and 
  error estimates. Scaling by powers of the radix does not cause rounding errors unless the result
  underflows or overflows. Rounding errors during scaling lead to refining with a matrix that is 
  not equivalent to the input matrix, producing error estimates that may not be reliable.

B:  array, dimension (LDB,NRHS). On entry, the N-by-NRHS right hand side matrix B. On exit,
  if EQUED = 'N', B is not modified; if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by 
  diag(R)*B; if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is overwritten by diag(C)*B.

LDB: The leading dimension of the array B.  LDB >= max(1,N).

X:  array, dimension (LDX,NRHS). If INFO = 0, the N-by-NRHS solution matrix X to the original 
  system of equations.  Note that A and B are modified on exit if EQUED .ne. 'N', and the 
  solution to the equilibrated system is   
    inv(diag(C))*X if TRANS = 'N' and EQUED = 'C' or 'B', or
    inv(diag(R))*X if TRANS = 'T' or 'C' and EQUED = 'R' or 'B'.

LDX:  The leading dimension of the array X.  LDX >= max(1,N).

RCOND:  Reciprocal scaled condition number.  This is an estimate of the reciprocal Skeel condition 
  number of the matrix A after equilibration (if done).  If this is less than the machine precision 
  (in particular, if it is zero), the matrix is singular to working precision. Note that the error 
  may still be small even if this number is very small and the matrix appears ill-conditioned.

RPVGRW:  Reciprocal pivot growth.  On exit, this contains the reciprocal pivot growth factor 
  norm(A)/norm(U). The "max absolute element" norm is used.  If this is much less than 1, then the
  stability of the LU factorization of the (equilibrated) matrix A could be poor. This also means 
  that the solution X, estimated condition numbers, and error bounds could be unreliable. If 
  factorization fails with 0<INFO<=N, then this contains the reciprocal pivot growth factor for the
  leading INFO columns of A.  In DGESVX, this quantity is returned in WORK(1).

BERR: array, dimension (NRHS). Componentwise relative backward error. This is the componentwise 
  relative backward error of each solution vector X(j) (i.e., the smallest relative change in any 
  element of A or B that makes X(j) an exact solution).

N_ERR_BNDS: Number of error bounds to return for each right hand side and each type (normwise or 
  componentwise).  See ERR_BNDS_NORM and ERR_BNDS_COMP below.

ERR_BNDS_NORM: array, dimension (NRHS, N_ERR_BNDS). For each right-hand side, this array contains
  information about various error bounds and condition numbers corresponding to the normwise 
  relative error, which is defined as follows: Normwise relative error in the ith solution vector:

             max_j (abs(XTRUE(j,i) - X(j,i)))
            ------------------------------
                  max_j abs(X(j,i))

  The array is indexed by the type of error information as described below. There currently are up 
  to three pieces of information returned. The first index in ERR_BNDS_NORM(i,:) corresponds to the
  ith right-hand side. The second index in ERR_BNDS_NORM(:,err) contains the following three 
  fields:
    err = 1 "Trust/don't trust" boolean. Trust the answer if the reciprocal condition number is 
             less than the threshold sqrt(n) * dlamch('Epsilon').
    err = 2 "Guaranteed" error bound: The estimated forward error, almost certainly within a factor
            of 10 of the true error so long as the next entry is greater than the threshold 
            sqrt(n) * dlamch('Epsilon'). This error bound should only be trusted if the previous 
            boolean is true.
    err = 3 Reciprocal condition number: Estimated normwise reciprocal condition number. Compared
            with the threshold sqrt(n) * dlamch('Epsilon') to determine if the error estimate is 
            "guaranteed". These reciprocal condition numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) 
            for some appropriately scaled matrix Z. Let Z = S*A, where S scales each row by a 
            power of the radix so all absolute row sums of Z are approximately 1.
    See Lapack Working Note 165 for further details and extra  cautions.

ERR_BNDS_COMP: array, dimension (NRHS, N_ERR_BNDS). For each right-hand side, this array contains 
  information about various error bounds and condition numbers corresponding to the componentwise 
  relative error, which is defined as follows: Componentwise relative error in the ith solution 
  vector:
                    abs(XTRUE(j,i) - X(j,i))
             max_j ----------------------
                         abs(X(j,i))

  The array is indexed by the right-hand side i (on which the componentwise relative error 
  depends), and the type of error information as described below. There currently are up to three
  pieces of information returned for each right-hand side. If componentwise accuracy is not 
  requested (PARAMS(3) = 0.0), then ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at 
  most the first (:,N_ERR_BNDS) entries are returned. The first index in ERR_BNDS_COMP(i,:) 
  corresponds to the ith right-hand side.
  The second index in ERR_BNDS_COMP(:,err) contains the following three fields:
   err = 1 "Trust/don't trust" boolean. Trust the answer if the reciprocal condition number is 
           less than the threshold sqrt(n) * dlamch('Epsilon').
   err = 2 "Guaranteed" error bound: The estimated forward error, almost certainly within a factor
           of 10 of the true error so long as the next entry is greater than the threshold
           sqrt(n) * dlamch('Epsilon'). This error bound should only be trusted if the previous 
           boolean is true.
   err = 3 Reciprocal condition number: Estimated componentwise reciprocal condition number. 
           Compared with the threshold sqrt(n) * dlamch('Epsilon') to determine if the error 
           estimate is "guaranteed". These reciprocal condition 
           numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some appropriately scaled matrix Z.
           Let Z = S*(A*diag(x)), where x is the solution for the current right-hand side and S 
           scales each row of A*diag(x) by a power of the radix so all absolute row sums of Z are 
           approximately 1.
   See Lapack Working Note 165 for further details and extra cautions.


NPARAMS:  Specifies the number of parameters set in PARAMS.  If .LE. 0, the PARAMS array is never 
  referenced and default values are used.

PARAMS:  array, dimension (NPARAMS). Specifies algorithm parameters.  If an entry is .LT. 0.0, then
  that entry will be filled with default value used for that parameter.  Only positions up to 
  NPARAMS are accessed; defaults are used for higher-numbered parameters.
  PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative refinement or not.
    Default: 1.0D+0
    = 0.0 : No refinement is performed, and no error bounds are computed.
    = 1.0 : Use the extra-precise refinement algorithm. (other values are reserved for future use)
  PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual computations allowed for refinement.
    Default: 10
    Aggressive: Set to 100 to permit convergence using approximate factorizations or 
    factorizations other than LU. If the factorization uses a technique other than Gaussian 
    elimination, the guarantees in err_bnds_norm and err_bnds_comp may no longer be trustworthy.
  PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code will attempt to find a solution with
    small componentwise relative error in the double-precision algorithm. Positive is true, 0.0 is 
    false. Default: 1.0 (attempt componentwise convergence)

WORK:  array, dimension (4*N)
IWORK: dimension (N)
INFO: 
  = 0:  Successful exit. The solution to every right-hand side is guaranteed.
  < 0:  If INFO = -i, the i-th argument had an illegal value
  > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization has been completed, but the 
      factor U is exactly singular, so the solution and error bounds could not be computed. 
      RCOND = 0 is returned.
  = N+J: The solution corresponding to the Jth right-hand side is not guaranteed. The solutions 
    corresponding to other right-hand sides K with K > J may not be guaranteed as well, but only 
    the first such right-hand side is reported. If a small componentwise error is not requested 
    (PARAMS(3) = 0.0) then the Jth right-hand side is the first with a normwise error bound that is
    not guaranteed (the smallest J such that ERR_BNDS_NORM(J,1) = 0.0). By default 
    (PARAMS(3) = 1.0) the Jth right-hand side is the first with either a normwise or componentwise 
    error bound that is not guaranteed (the smallest  J such that either ERR_BNDS_NORM(J,1) = 0.0 
    or ERR_BNDS_COMP(J,1) = 0.0). See the definition of ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). 
    To get information about all of the right-hand sides check ERR_BNDS_NORM or ERR_BNDS_COMP. */
template<class T>
int gbsvxx(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab, 
  integer *ldab, T *afb, integer *ldafb, integer *ipiv, char *equed, T *r__, T *c__, T *b, 
  integer *ldb, T *x, integer *ldx, T *rcond, T *rpvgrw, T *berr, integer *n_err_bnds__, 
  T *err_bnds_norm__, T *err_bnds_comp__, integer *nparams, T *params, T *work, integer *iwork, 
  integer *info, ftnlen fact_len, ftnlen trans_len, ftnlen equed_len);


//=================================================================================================
// Computational routines:


//-------------------------------------------------------------------------------------------------

/**
Purpose:
DGBCON estimates the reciprocal of the condition number of a real general band matrix A, in either 
the 1-norm or the infinity-norm, using the LU factorization computed by DGBTRF. An estimate is 
obtained for norm(inv(A)), and the reciprocal of the condition number is computed as
RCOND = 1 / ( norm(A) * norm(inv(A)) ).

Arguments:
NORM:  Specifies whether the 1-norm condition number or the infinity-norm condition number is required:
       = '1' or 'O':  1-norm;
       = 'I':         Infinity-norm.
N:     The order of the matrix A.  N >= 0.
KL:    The number of subdiagonals within the band of A.  KL >= 0.
KU:    The number of superdiagonals within the band of A.  KU >= 0.
AB:    array, dimension (LDAB,N). Details of the LU factorization of the band matrix A, as computed
       by DGBTRF.  U is stored as an upper triangular band matrix with KL+KU superdiagonals in 
       rows 1 to KL+KU+1, and the multipliers used during the factorization are stored in 
       rows KL+KU+2 to 2*KL+KU+1.
LDAB:  The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
IPIV:  array, dimension (N). The pivot indices; for 1 <= i <= N, row i of the matrix was 
       interchanged with row IPIV(i).
ANORM: If NORM = '1' or 'O', the 1-norm of the original matrix A.
       If NORM = 'I', the infinity-norm of the original matrix A.
RCOND: The reciprocal of the condition number of the matrix A, computed as 
       RCOND = 1/(norm(A) * norm(inv(A))).
WORK:  array, dimension (3*N)
IWORK: array, dimension (N).
INFO:  = 0:  successful exit
       < 0: if INFO = -i, the i-th argument had an illegal value  */
template<class T>
int gbcon(char *norm, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, integer *ipiv, 
  T *anorm, T *rcond, T *work, integer *iwork, integer *info, ftnlen norm_len);



//-------------------------------------------------------------------------------------------------

/**

Purpose:
DGBEQU computes row and column scalings intended to equilibrate an M-by-N band matrix A and reduce
its condition number.  R returns the row scale factors and C the column scale factors, chosen to 
try to make the largest element in each row and column of the matrix B with 
elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1. R(i) and C(j) are restricted to be between
SMLNUM = smallest safe number and BIGNUM = largest safe number.  Use of these scaling factors is 
not guaranteed to reduce the condition number of A but works well in practice.

Arguments:
M:      The number of rows of the matrix A.  M >= 0.
N:      The number of columns of the matrix A.  N >= 0.
KL:     The number of subdiagonals within the band of A.  KL >= 0.
KU:     The number of superdiagonals within the band of A.  KU >= 0.
AB:     array, dimension (LDAB,N). The band matrix A, stored in rows 1 to KL+KU+1. The j-th column
        of A is stored in the j-th column of the array AB as follows:
          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl).
LDAB:   The leading dimension of the array AB.  LDAB >= KL+KU+1.
R:      array, dimension (M). If INFO = 0, or INFO > M, R contains the row scale factors for A.
C:      array, dimension (N). If INFO = 0, C contains the column scale factors for A.

ROWCND: If INFO = 0 or INFO > M, ROWCND contains the ratio of the smallest R(i) to the largest 
        R(i).  If ROWCND >= 0.1 and AMAX is neither too large nor too small, it is not worth 
        scaling by R.
COLCND: If INFO = 0, COLCND contains the ratio of the smallest C(i) to the largest C(i).  If 
        COLCND >= 0.1, it is not worth scaling by C.
AMAX:   Absolute value of largest matrix element.  If AMAX is very close to overflow or very close
        to underflow, the matrix should be scaled.
INFO:   = 0:  successful exit
        < 0:  if INFO = -i, the i-th argument had an illegal value
        > 0:  if INFO = i, and i is
        <= M:  the i-th row of A is exactly zero
        >  M:  the (i-M)-th column of A is exactly zero  */
template<class T>
int gbequ(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *r__, T *c__, 
  T *rowcnd, T *colcnd, T *amax, integer *info);

//-------------------------------------------------------------------------------------------------

/**
Purpose:
DGBRFS improves the computed solution to a system of linear equations when the coefficient matrix 
is banded, and provides error bounds and backward error estimates for the solution.

Arguments:
TRANS: Specifies the form of the system of equations:
       = 'N':  A * X = B     (No transpose)
       = 'T':  A**T * X = B  (Transpose)
       = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
N:     The order of the matrix A.  N >= 0.
KL:    The number of subdiagonals within the band of A.  KL >= 0.
KU:    The number of superdiagonals within the band of A.  KU >= 0.
NRHS:  The number of right hand sides, i.e., the number of columns of the matrices B and X.  
       NRHS >= 0.
AB:    array, dimension (LDAB,N). The original band matrix A, stored in rows 1 to KL+KU+1. The j-th
       column of A is stored in the j-th column of the array AB as follows:
       AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).
LDAB:  The leading dimension of the array AB.  LDAB >= KL+KU+1.
AFB:   array, dimension (LDAFB,N). Details of the LU factorization of the band matrix A, as 
       computed by DGBTRF.  U is stored as an upper triangular band matrix with KL+KU 
       superdiagonals in rows 1 to KL+KU+1, and the multipliers used during the factorization are 
       stored in rows KL+KU+2 to 2*KL+KU+1.
LDAFB: The leading dimension of the array AFB.  LDAFB >= 2*KL*KU+1.
IPIV:  array, dimension (N). The pivot indices from DGBTRF; for 1<=i<=N, row i of the matrix was 
       interchanged with row IPIV(i).
B:     array, dimension (LDB,NRHS). The right hand side matrix B.
LDB:   The leading dimension of the array B.  LDB >= max(1,N).
X:     array, dimension (LDX,NRHS). On entry, the solution matrix X, as computed by DGBTRS. On 
       exit, the improved solution matrix X.
LDX:   The leading dimension of the array X.  LDX >= max(1,N).
FERR:  array, dimension (NRHS). The estimated forward error bound for each solution vector X(j) 
       (the j-th column of the solution matrix X). If XTRUE is the true solution corresponding to 
       X(j), FERR(j) is an estimated upper bound for the magnitude of the largest element in 
       (X(j) - XTRUE) divided by the magnitude of the largest element in X(j).  The estimate is as 
       reliable as the estimate for RCOND, and is almost always a slight overestimate of the true 
       error.
BERR:  array, dimension (NRHS). The componentwise relative backward error of each solution vector 
       X(j) (i.e., the smallest relative change in any element of A or B that makes X(j) an exact 
       solution).
WORK:  array, dimension (3*N)
IWORK: dimension (N)
INFO:  = 0:  successful exit
       < 0:  if INFO = -i, the i-th argument had an illegal value

Internal Parameters:
ITMAX: is the maximum number of steps of iterative refinement. */
template<class T>
int gbrfs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab, integer *ldab, 
  T *afb, integer *ldafb, integer *ipiv, T *b, integer *ldb, T *x, integer *ldx, T *ferr, T *berr, 
  T *work, integer *iwork, integer *info, ftnlen trans_len);

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


/** Returns true, if x is NaN, flase otherwise (re-implemented by Robin Schmidt). */
template<class T>
logical isnan(T *x) { return *x != *x; }

//-------------------------------------------------------------------------------------------------

/**
Purpose:
DLACN2 estimates the 1-norm of a square, real matrix A. Reverse communication is used for 
evaluating matrix-vector products.

Arguments:
N:     The order of the matrix.  N >= 1.
V:     array, dimension (N). On the final return, V = A*W,  where  EST = norm(V)/norm(W) 
       (W is not returned).
X:     array, dimension (N). On an intermediate return, X should be overwritten by 
       A    * X,  if KASE=1,
       A**T * X,  if KASE=2,
       and DLACN2 must be re-called with all the other parameters unchanged.
ISGN:  array, dimension (N)
EST:   On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be unchanged from the previous call
       to DLACN2. On exit, EST is an estimate (a lower bound) for norm(A).
KASE:  On the initial call to DLACN2, KASE should be 0. On an intermediate return, KASE will be 
       1 or 2, indicating whether X should be overwritten by A * X  or A**T * X. On the final return
       from DLACN2, KASE will again be 0.        
ISAVE: is INTEGER array, dimension (3). ISAVE is used to save variables between calls to DLACN2

Further Details:
Originally named SONEST, dated March 16, 1988. This is a thread safe version of DLACON, which uses
the array ISAVE in place of a SAVE statement, as follows:
   DLACON     DLACN2
    JUMP     ISAVE(1)
    J        ISAVE(2)
    ITER     ISAVE(3)  */
template<class T>
int lacn2(integer *n, T *v, T *x, integer *isgn, T *est, integer *kase, integer *isave);

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


/** DLAMCH determines double precision machine parameters.

Arguments:
CMACH: Specifies the value to be returned by DLAMCH:
       = 'E' or 'e',   DLAMCH := eps       = relative machine precision
       = 'S' or 's ,   DLAMCH := sfmin     = safe minimum, such that 1/sfmin does not overflow
       = 'B' or 'b',   DLAMCH := base      = base of the machine
       = 'P' or 'p',   DLAMCH := eps*base  = eps*base
       = 'N' or 'n',   DLAMCH := t         = number of (base) digits in the mantissa
       = 'R' or 'r',   DLAMCH := rnd       = 1.0 when rounding occurs in addition, 0.0 otherwise
       = 'M' or 'm',   DLAMCH := emin      = minimum exponent before (gradual) underflow
       = 'U' or 'u',   DLAMCH := rmin      = underflow threshold - base**(emin-1)
       = 'L' or 'l',   DLAMCH := emax      = largest exponent before overflow
       = 'O' or 'o',   DLAMCH := rmax      = overflow threshold  - (base**emax)*(1-eps) */
doublereal lamch(char *cmach, ftnlen cmach_len);
//template<class T> ..hmm - can't be templatized bcs the datatype appears only in the return value
// maybe pass a dummy variable - but that requires all calling code to be adapted - but maybe that
// function is only called in a few palces, so that doesn't matter? we'll see. -if it's practical,
// get rid of the "d" - it should be only lamch 
// i think,...we need to somehow call std::numeric_limits<dummy>::eps, etc. in the implementation 
// to make it work ...or use function templates machineEps, machineSafeMin and provide explicit
// specializations or better: numericEps,

// For this routine, i need to actually to actually (re)write some code - implementthe functions
// that return the requested values for the machine precision - using explicit intantiations

//-------------------------------------------------------------------------------------------------

/**
Purpose:
DLANTB  returns the value of the one norm,  or the Frobenius norm, or the  infinity norm, or the 
element of  largest absolute value  of an n by n triangular band matrix A,  with ( k + 1 ) diagonals.
  DLANTB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
           (
           ( norm1(A),         NORM = '1', 'O' or 'o'
           (
           ( normI(A),         NORM = 'I' or 'i'
           (
           ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
where  norm1  denotes the  one norm of a matrix (maximum column sum), normI  denotes the  infinity 
norm of a matrix  (maximum row sum) and normF  denotes the  Frobenius norm of a matrix (square root
of sum of squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.

Arguments:
NORM:  Specifies the value to be returned in DLANTB as described above.
UPLO:  Specifies whether the matrix A is upper or lower triangular.
       = 'U':  Upper triangular
       = 'L':  Lower triangular
DIAG:  Specifies whether or not the matrix A is unit triangular.
       = 'N':  Non-unit triangular
       = 'U':  Unit triangular
N:     The order of the matrix A.  N >= 0.  When N = 0, DLANTB is set to zero.
K:     The number of super-diagonals of the matrix A if UPLO = 'U', or the number of sub-diagonals 
       of the matrix A if UPLO = 'L'. K >= 0.
AB:    array, dimension (LDAB,N). The upper or lower triangular band matrix A, stored in the first 
       k+1 rows of AB.  The j-th column of A is stored in the j-th column of the array AB as follows:
       if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j; 
       if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).
       Note that when DIAG = 'U', the elements of the array AB corresponding to the diagonal 
       elements of the matrix A are not referenced, but are assumed to be one.
LDAB:  The leading dimension of the array AB.  LDAB >= K+1.
WORK:  array, dimension (MAX(1,LWORK)), where LWORK >= N when NORM = 'I'; otherwise, WORK is not
       referenced.  */
template<class T>
T lantb(char *norm, char *uplo, char *diag, integer *n, integer *k, T *ab, integer *ldab, T *work, 
  ftnlen norm_len, ftnlen uplo_len, ftnlen diag_len);

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

/**
Purpose:
DLASSQ  returns the values  scl  and  smsq  such that

  ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,

where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is assumed to be non-negative and
scl  returns the value 

   scl = max( scale, abs( x( i ) ) ).

scale and sumsq must be supplied in SCALE and SUMSQ and scl and smsq are overwritten on SCALE 
and SUMSQ respectively. he routine makes only one pass through the vector x.

Arguments:
N:     The number of elements to be used from the vector X.
X:     array, dimension (N). The vector for which a scaled sum of squares is computed.
       x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
INCX:  The increment between successive values of the vector X. INCX > 0.
SCALE: On entry, the value  scale  in the equation above. On exit, SCALE is overwritten with scl, 
       the scaling factor for the sum of squares.
SUMSQ: On entry, the value  sumsq  in the equation above. On exit, SUMSQ is overwritten with smsq, 
       the basic sum of squares from which  scl  has been factored out. */
template<class T>
int lassq(integer *n, T *x, integer *incx, T *scale, T *sumsq);


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

//-------------------------------------------------------------------------------------------------

/**
Purpose:
DLATBS solves one of the triangular systems

  A *x = s*b  or  A**T*x = s*b

with scaling to prevent overflow, where A is an upper or lower triangular band matrix. Here A**T
denotes the transpose of A, x and b are n-element vectors, and s is a scaling factor, usually less 
than or equal to 1, chosen so that the components of x will be less than the overflow threshold. 
If the unscaled problem will not cause overflow, the Level 2 BLAS routine DTBSV is called. If the
matrix A is singular (A(j,j) = 0 for some j), then s is set to 0 and a non-trivial solution to 
A*x = 0 is returned.

Arguments:
UPLO:   Specifies whether the matrix A is upper or lower triangular.
        = 'U':  Upper triangular
        = 'L':  Lower triangular
TRANS:  Specifies the operation applied to A.
        = 'N':  Solve A * x = s*b  (No transpose)
        = 'T':  Solve A**T* x = s*b  (Transpose)
        = 'C':  Solve A**T* x = s*b  (Conjugate transpose = Transpose)
DIAG:   Specifies whether or not the matrix A is unit triangular.
        = 'N':  Non-unit triangular
        = 'U':  Unit triangular
NORMIN: Specifies whether CNORM has been set or not.
        = 'Y':  CNORM contains the column norms on entry
        = 'N':  CNORM is not set on entry. On exit, the norms will be computed and stored in CNORM.
N:      The order of the matrix A.  N >= 0.
KD:     The number of subdiagonals or superdiagonals in the triangular matrix A.  KD >= 0.
AB:     array, dimension (LDAB,N). The upper or lower triangular band matrix A, stored in the first
        KD+1 rows of the array. The j-th column of A is stored in the j-th column of the array AB 
        as follows:
        if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
        if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
LDAB:   The leading dimension of the array AB.  LDAB >= KD+1.
X:      array, dimension (N). On entry, the right hand side b of the triangular system. On exit, X 
        is overwritten by the solution vector x.
SCALE:  The scaling factor s for the triangular system
             A * x = s*b  or  A**T* x = s*b.
        If SCALE = 0, the matrix A is singular or badly scaled, and the vector x is an exact or 
        approximate solution to A*x = 0.
CNORM:  array, dimension (N). If NORMIN = 'Y', CNORM is an input argument and CNORM(j) contains the 
        norm of the off-diagonal part of the j-th column of A. If TRANS = 'N', CNORM(j) must be 
        greater than or equal to the infinity-norm, and if TRANS = 'T' or 'C', CNORM(j) must be 
        greater than or equal to the 1-norm. If NORMIN = 'N', CNORM is an output argument and 
        CNORM(j) returns the 1-norm of the offdiagonal part of the j-th column of A.
INFO:   = 0:  successful exit
        < 0:  if INFO = -k, the k-th argument had an illegal value

Further Details:
A rough bound on x is computed; if that is less than overflow, DTBSV is called, otherwise, specific
code is used which checks for possible overflow or divide-by-zero at every operation. A columnwise 
scheme is used for solving A*x = b. The basic algorithm if A is lower triangular is

  x[1:n] := b[1:n]
  for j = 1, ..., n
       x(j) := x(j) / A(j,j)
       x[j+1:n] := x[j+1:n] - x(j) * A[j+1:n,j]
  end

Define bounds on the components of x after j iterations of the loop:

  M(j) = bound on x[1:j]
  G(j) = bound on x[j+1:n]

Initially, let M(0) = 0 and G(0) = max{x(i), i=1,...,n}. Then for iteration j+1 we have

  M(j+1) <= G(j) / | A(j+1,j+1) |
  G(j+1) <= G(j) + M(j+1) * | A[j+2:n,j+1] |
         <= G(j) ( 1 + CNORM(j+1) / | A(j+1,j+1) | )

where CNORM(j+1) is greater than or equal to the infinity-norm of column j+1 of A, not counting the 
diagonal. Hence

  G(j) <= G(0) product ( 1 + CNORM(i) / | A(i,i) | )
               1<=i<=j

 |x(j)| <= ( G(0) / |A(j,j)| ) product ( 1 + CNORM(i) / |A(i,i)| )
                               1<=i< j

Since |x(j)| <= M(j), we use the Level 2 BLAS routine DTBSV if the reciprocal of the largest 
M(j), j=1,..,n, is larger than max(underflow, 1/overflow). The bound on x(j) is also used to 
determine when a step in the columnwise method can be performed without fear of overflow. If the 
computed bound is greater than a large constant, x is scaled to prevent overflow, but if the bound
overflows, x is set to 0, x(j) to 1, and scale to 0, and a non-trivial solution to A*x = 0 is 
found. Similarly, a row-wise scheme is used to solve A**T*x = b. The basic algorithm for A upper 
triangular is

  for j = 1, ..., n
       x(j) := ( b(j) - A[1:j-1,j]**T * x[1:j-1] ) / A(j,j)
  end

We simultaneously compute two bounds
  G(j) = bound on ( b(i) - A[1:i-1,i]**T * x[1:i-1] ), 1<=i<=j
  M(j) = bound on x(i), 1<=i<=j
The initial values are G(0) = 0, M(0) = max{b(i), i=1,..,n}, and we add the constraint 
G(j) >= G(j-1) and M(j) >= M(j-1) for j >= 1. Then the bound on x(j) is

  M(j) <= M(j-1) * ( 1 + CNORM(j) ) / | A(j,j) |
       <= M(0) * product ( ( 1 + CNORM(i) ) / |A(i,i)| )
                 1<=i<=j

and we can safely call DTBSV if 1/M(n) and 1/G(n) are both greater than 
max(underflow, 1/overflow). */
template<class T>
int latbs(char *uplo, char *trans, char *diag, char *normin, integer *n, integer *kd, T *ab, 
  integer *ldab, T *x, T *scale, T *cnorm, integer *info, ftnlen uplo_len, ftnlen trans_len, 
  ftnlen diag_len, ftnlen normin_len);







// maybe the LaPack source files should be split into LaPackGB, etc. or: LaPackDrivers, 
// LapackComputation, LaPackAuxiliary


}