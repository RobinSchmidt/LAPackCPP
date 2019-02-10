#pragma once

#include "../GeneratedByF2C/f2c.h" // todo: move file to this directory

namespace BlasCPP {


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



//=================================================================================================

/** \name BLAS level 2 routines (operations involving matrices and vectors) */


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
INCY:  On entry, INCY specifies the increment for the elements of Y. INCY must not be zero. */
template<class T>
int gbmv(char *trans, integer *m, integer *n, integer *kl, integer *ku, T *alpha, T *a, 
  integer *lda, T *x, integer *incx, T *beta, T *y, integer *incy, ftnlen trans_len);



}