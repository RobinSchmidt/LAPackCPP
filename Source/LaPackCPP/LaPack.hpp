#pragma once

// maybe we will later need to include Blas.hpp
#include "../GeneratedByF2C/f2c.h" // todo: move file to this directory

namespace LaPackCPP {

//-------------------------------------------------------------------------------------------------

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



// maybe the LaPack source files should be split into LaPackGB, etc.


}