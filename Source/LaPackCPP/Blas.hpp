#pragma once

#include "../GeneratedByF2C/f2c.h" // todo: move file to libf2c folder

namespace BlasCPP {


// How can we make sure, that client code can use a drop-in replacement for these (unoptimized)
// reference implementations of the BLAS routines? Not allowing that would go totally against the
// design idea of LAPACK (which is to rely on lower-level linear algebra routines that can be
// optimized for particular target platform and then dropped in). Maybe the LAPACK routines call
// a general template function like axpy (not daxpy, saxpy, caxpy or zaxpy)...and then it must be
// somehow dispatched at compile-time which of the 4 (d,s,c,z) is chosen depending on the template
// parameter for axpy. The design goal of this C++ translation is to NOT have separate routines for
// different datatypes and instead just have one axpy function...maybe for some things that are 
// specific to complex numbers, that may not work out completely - we'll see - but for many of the 
// BLAS and LAPACK routines, it should be possible to write them in a type independent way as 
// templates. Maybe an explicit specialization of axpy for double should call daxpy and that daxpy
// may then be defined somewhere else. If we don't want to use such an explicit specialization, we 
// just request an explicit instantiation of axpy from the compiler - this can be done by letting
// the user decide which source files are compiled with the application - whether it's the ones 
// with explicit instantiations or the ones with explicit specializations...something along those
// lines...
// Maybe i should implement my own, very naive, BLAS routines (no loop unrolling, etc.). That stuff
// might better be left to the compiler anyway - maybe unrolling by 4 (as this reference 
// implementation does) is not optimal and an optimizing compiler would unroll by some other number 
// (8,16,..) instead? So maybe a more naive BLAS implementation, in addition to be more readable, 
// could indeed perfom better? -> try it! With such a "NaiveBlas" in place, we could also try the 
// replacement mechanism, once it's implemented, and swap between the reference implementation 
// maybe "RefBlas" and make some comparison benchmarks. ..also try https://www.netlib.org/atlas/

//=================================================================================================

/** \name BLAS level 1 routines (operations involving scalars and vectors) */

/** Computes constant times a vector plus a vector: y = a*x + y
N:    Number of elements in input vector(s)
a:    On entry, a specifies the scalar alpha.
x:    array, dimension ( 1 + ( N - 1 )*abs( incX ) ), input
incX: storage spacing between elements of x
y:    array, dimension ( 1 + ( N - 1 )*abs( incY ) ), input/output
incY: storage spacing between elements of y  */
template<class T>
int axpy(long int* N, T* a, T* x, long int* incX, T* y, long int* incY);

/** LSAME returns .TRUE. if CA is the same letter as CB regardless of case. CA and CB specify the 
single characters to be compared. */
logical lsame(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len);



/** xerbla is an error handler

Purpose:
XERBLA is an error handler for the LAPACK routines. It is called by an LAPACK routine if an input 
parameter has an invalid value. A message is printed and execution stops. Installers may consider 
modifying the STOP statement in order to call system-specific exception-handling facilities. 

Arguments:
SRNAME: The name of the routine which called XERBLA.
INFO:   The position of the invalid parameter in the parameter list of the calling routine. */
int xerbla(char *srname, integer *info, ftnlen srname_len);


//=================================================================================================

/** \name BLAS level 2 routines (operations involving matrices and vectors) */

/** ger updates matrix A := alpha*x*y**T + A

Purpose: 
ger performs the rank 1 operation 
A := alpha*x*y**T + A
where alpha is a scalar, x is an m element vector, y is an n element vector and A is an m by n 
matrix.

Arguments:
M:     On entry, M specifies the number of rows of the matrix A. M must be at least zero.
N:     On entry, N specifies the number of columns of the matrix A. N must be at least zero.
ALPHA: On entry, ALPHA specifies the scalar alpha.
X:     array, dimension at least ( 1 + ( m - 1 )*abs( INCX ) ). Before entry, the incremented array
       X must contain the m element vector x.
INCX:  On entry, INCX specifies the increment for the elements of X. INCX must not be zero.
Y:     array, dimension at least ( 1 + ( n - 1 )*abs( INCY ) ). Before entry, the incremented array 
       Y must contain the n element vector y.
INCY:  On entry, INCY specifies the increment for the elements of Y. INCY must not be zero.
A:     array, dimension ( LDA, N ). Before entry, the leading m by n part of the array A must 
       contain the matrix of coefficients. On exit, A is overwritten by the updated matrix.
LDA:   On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. 
       LDA must be at least max( 1, m ). */
int ger(integer *m, integer *n, doublereal *alpha, doublereal *x, integer *incx, doublereal *y, 
  integer *incy, doublereal *a, integer *lda);
// templatize!

//-------------------------------------------------------------------------------------------------

/** gbmv performs general banded matrix-vector multiplication

Purpose:
gbmv performs one of the matrix-vector operations
y := alpha*A*x    + beta*y,   or
y := alpha*A**T*x + beta*y
where alpha and beta are scalars, x and y are vectors and A is an m by n band matrix, with kl 
sub-diagonals and ku super-diagonals.

Arguments: 
TRANS: On entry, TRANS specifies the operation to be performed as follows: 
       TRANS = 'N' or 'n'   y := alpha*A*x + beta*y
       TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y. 
       TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y. 
M:     On entry, M specifies the number of rows of the matrix A. M must be at least zero. 
N:     On entry, N specifies the number of columns of the matrix A. N must be at least zero. 
KL:    On entry, KL specifies the number of sub-diagonals of the matrix A. KL must satisfy 
       0 .le. KL. 
KU:    On entry, KU specifies the number of super-diagonals of the matrix A. KU must satisfy  
       0 .le. KU. 
ALPHA: On entry, ALPHA specifies the scalar alpha. 
A:     array, dimension ( LDA, N ). Before entry, the leading ( kl + ku + 1 ) by n part of the 
       array A must contain the matrix of coefficients, supplied column by column, with the leading
       diagonal of the matrix in row ( ku + 1 ) of the array, the first super-diagonal starting at 
       position 2 in row ku, the first sub-diagonal starting at position 1 in row ( ku + 2 ), and 
       so on. Elements in the array A that do not correspond to elements in the band matrix (such 
       as the top left ku by ku triangle) are not referenced. The following program segment will 
       transfer a band matrix from conventional full matrix storage to band storage: 
          DO 20, J = 1, N 
             K = KU + 1 - J 
             DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL ) 
                A( K + I, J ) = matrix( I, J ) 
       10    CONTINUE 
       20 CONTINUE 
LDA:   On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. LDA 
       must be at least ( kl + ku + 1 ). 
X:     array, dimension at least ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' and at least
       ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. Before entry, the incremented array X must contain the 
       vector x. 
INCX:  On entry, INCX specifies the increment for the elements of X. INCX must not be zero. 
BETA:  On entry, BETA specifies the scalar beta. When BETA is supplied as zero then Y need not be 
       set on input. 
Y:     array, dimension at least ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' and at least 
       ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. Before entry, the incremented array Y must contain the 
       vector y. On exit, Y is overwritten by the updated vector y.   
INCY:  On entry, INCY specifies the increment for the elements of Y. INCY must not be zero. 

Robin's notes: 
trans_len is not used in the routine and the returned value is always zero (so, it has no meaning). */
template<class T>
int gbmv(char *trans, integer *m, integer *n, integer *kl, integer *ku, T *alpha, T *a, 
  integer *lda, T *x, integer *incx, T *beta, T *y, integer *incy, ftnlen trans_len);


//=================================================================================================

/** \name BLAS level 3 routines (operations involving two matrices) */

/**

Purpose:
GEMM  performs one of the matrix-matrix operations
 C := alpha*op( A )*op( B ) + beta*C
where  op( X ) is one of
 op( X ) = X   or   op( X ) = X**T,
alpha and beta are scalars, and A, B and C are matrices, with op( A ) an m by k matrix, op( B )  
a k by n matrix and  C an m by n matrix. 

Arguments:
TRANSA: On entry, TRANSA specifies the form of op( A ) to be used in the matrix multiplication as 
        follows: 
        TRANSA = 'N' or 'n',  op( A ) = A. 
        TRANSA = 'T' or 't',  op( A ) = A**T. 
        TRANSA = 'C' or 'c',  op( A ) = A**T. 
TRANSB: On entry, TRANSB specifies the form of op( B ) to be used in the matrix multiplication as 
        follows
        TRANSB = 'N' or 'n',  op( B ) = B.
        TRANSB = 'T' or 't',  op( B ) = B**T.
        TRANSB = 'C' or 'c',  op( B ) = B**T.
M:      On entry,  M  specifies  the number  of rows  of the  matrix op( A )  and of the  matrix C.  
        M  must  be at least  zero.
N:      On entry,  N  specifies the number  of columns of the matrix op( B ) and the number of 
        columns of the matrix C. N must be at least zero.
K:      On entry,  K  specifies  the number of columns of the matrix op( A ) and the number of rows 
        of the matrix op( B ). K must be at least  zero.
ALPHA:  On entry, ALPHA specifies the scalar alpha.
A:      array, dimension ( LDA, ka ), where ka is k  when  TRANSA = 'N' or 'n',  and is m otherwise.
        Before entry with  TRANSA = 'N' or 'n',  the leading  m by k part of the array A must contain 
        the matrix A, otherwise the leading  k by m  part of the array  A  must contain  the matrix A.
LDA:    On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. 
        When  TRANSA = 'N' or 'n' then LDA must be at least  max( 1, m ), otherwise  LDA must be at
        least  max( 1, k ).
B:      array, dimension ( LDB, kb ), where kb is n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
        Before entry with  TRANSB = 'N' or 'n',  the leading  k by n part of the array  B  must contain
        the matrix B, otherwise the leading n by k part of the array B must contain  the  matrix B.
LDB:    On entry, LDB specifies the first dimension of B as declared in the calling (sub) program. 
        When  TRANSB = 'N' or 'n' then LDB must be at least  max( 1, k ), otherwise  LDB must be at 
        least  max( 1, n ).
BETA:   On entry,  BETA  specifies the scalar  beta.  When  BETA  is supplied as zero then C need 
        not be set on input.
C:      array, dimension ( LDC, N ) Before entry, the leading  m by n  part of the array  C must 
        contain the matrix  C,  except when  beta  is zero, in which case C need not be set on 
        entry. On exit, the array  C  is overwritten by the  m by n  matrix 
        ( alpha*op( A )*op( B ) + beta*C ).
LDC:    On entry, LDC specifies the first dimension of C as declared in the calling (sub) program.   
        LDC must be at least max( 1, m ). */
int gemm(char *transa, char *transb, integer *m, integer *n, integer *k, doublereal *alpha, 
  doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *beta, doublereal *c__, 
  integer *ldc, ftnlen transa_len, ftnlen transb_len);
// templatize!


}