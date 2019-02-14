#pragma once

// maybe we will later need to include Blas.hpp
#include "../GeneratedByF2C/f2c.h" // todo: move file to libf2c folder

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