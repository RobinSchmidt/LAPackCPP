// LaPackTest1.cpp : Defines the entry point for the console application.
//
//#include "stdafx.h"

#include <vector>

#include "../../Source/Tests/TestUtils.cpp"

// todo: make a unity build and include file for all converted c code:
#include "../../Source/GeneratedByF2C/daxpy.c"

#include "../../Source/LaPackCPP/Blas.hpp"
using namespace BlasCPP;



std::vector<double> rangeVector(int N, double start, double inc)
{
  std::vector<double> v(N);
  for(size_t i = 0; i < N; i++)
    v[i] = start + i*inc;
  return v;
}

/** Utility class to wrap around an existing flat array of numbers to conveniently access matrix 
elements via row and column index via the () operator, taking row and column index as arguments. 
The access operator always takes as first index the row index and as second index the column index, 
regardless how the matrix is stored internally (column-major, in case of using this class). In math 
notation, the first index IS the row index and math notation does not care about storage formats 
and the access operator should just be consistent with math notation and abstract away from the 
storage format - that's its main purpose, actually */
template<class T>
class MatrixViewColumnMajor
{

  /** . */
  //T& operator()(int row, int column) {  }
  // todo: make a general MatrixView class that supports row-major and column-major element flat
  // access, access via row-pointers or column-pointers and 0-based and 1-based indexing (numerical 
  // recipies routines use 1-based (using pointer offset trickery) indexed, row-pointer access (at 
  // least mostly), LAPack routines use column-major, flat-array, 0-based access - maybe we can 
  // just do return (rowStride - base)*i + (columnStride - base)*j for generalized 
  // flat-array access? for row-major indexing, rowStride = numColumns and columnStride = 1, for
  // column major: rowStride = 1, columnStride = numRows, and base will be 0 for 0-based 
  // indexing and 1 for 1-based indexing - but we should make performance tests, if supporting all
  // these ways has an performance hit - the two multiplications can probably be done in parallel
  // so that may not matter - but the subtraction of 1 or 0 has to be in series, so it may have 
  // a performance impact - we'll see....- if it does make a difference, make different matrix 
  // classes for different access modes, so the access mode can be decided at compile time.
  // actually, it makes no sense anyway to be able to switch the access-mode at runtime - it makes 
  // sense only at compile time, for example, when we have translated an algorithm from a language 
  // with 1-based indexing and it is still in 1-based form - before eventually converting all 
  // accesses into 0-based form, it may be convenient to temporarily just be able to use 1-based
  // indexing


protected:

  int numRows, numColumns;
  T *data;

};


void mulMatVec(int M, int N, double* A[], double x[], double y[])
{

}



bool gbsvUnitTest()
{
  bool r = true;

  // We multiply the 9x9 band-diagonal matrix with 2 upper and 3 lower diagonals with 4 different 
  // vectors x (i.e. a 9x4 X-matrix) to obtain a B matrix. Then we solve the system A * X = B for
  // the solution matrix X and compare it to the original X matrix...

  // 11 12 13 00 00 00 00 00 00     1 9 1 0
  // 21 22 23 24 00 00 00 00 00     2 8 0 1
  // 31 32 33 34 35 00 00 00 00     3 7 2 0
  // 41 42 43 44 45 46 00 00 00     4 6 0 2
  // 00 52 53 54 55 56 57 00 00  *  5 5 3 0
  // 00 00 63 64 65 66 67 68 00     6 4 0 3
  // 00 00 00 74 75 76 77 78 79     7 3 4 0
  // 00 00 00 00 85 86 87 88 89     8 2 0 4
  // 00 00 00 00 00 96 97 98 99     9 1 5 0

  typedef const long int Int;  // for convenience
  static Int N    = 9;         // size of the system: NxN

  // the matrix itself for reference (lapack uses column-major indexing, so do we here):
  double refMat[N][N] =
  { { 11,21,31,41, 0, 0, 0, 0, 0 },    // 1st column (index 0)
    { 12,22,32,42,52, 0, 0, 0, 0 },    // 2nd column (index 1)
    { 13,23,33,43,53,63, 0, 0, 0 },    // etc.
    {  0,24,34,44,54,64,74, 0, 0 },
    {  0, 0,35,45,55,65,75,85, 0 },
    {  0, 0, 0,46,56,66,76,86,96 },
    {  0, 0, 0, 0,57,67,77,87,97 },
    {  0, 0, 0, 0, 0,68,78,88,98 },
    {  0, 0, 0, 0, 0, 0,79,89,99 } };  // 9th column (index 8)
  // maybe later store the matrix in LAPACK format to compare the general matrix multiply (gemv) 
  // result with the banded one (gbmv)



  // the banded storage format suitable for gbsv looks in this case like below, where column-major 
  // indexing is assumed (as always in LAPACK):

  // ** ** ** ** ** ++ ++ ++ ++
  // ** ** ** ** ++ ++ ++ ++ ++   fields with ++ are used internally by the algorithm, fields
  // ** ** ** ++ ++ ++ ++ ++ ++   with ** are not accessed
  // ** ** 13 24 35 46 57 68 79   upper diagonal 2
  // ** 12 23 34 45 56 67 78 89   upper diagonal 1
  // 11 22 33 44 55 66 77 88 99   main diagonal
  // 21 32 43 54 65 76 87 98 **   lower diagonal 1
  // 31 42 53 64 75 86 97 ** **   lower diagonal 2
  // 41 52 63 74 85 96 ** ** **   lower diagonal 3

  static Int kl   = 3;         // number of lower diagonals
  static Int ku   = 2;         // number of upper diagonals
  static Int nrhs = 4;         // number of right hand sides
  static Int ldab = 2*kl+ku+1; // leading dimension of ab >= 2*kl+ku+1 
  //static Int arab = 2*kl+ku+1; // allocated rows matrix ab
  static Int ldb  = N;         // leading dimension of b >= max(1,N). ...is N correct?
                               // rename to arb - allocated rows for b
  long int ipiv[N];            // pivot indices
  int info = 666;              // 0 on exit, if successful
  double _ = 0.0/sqrt(0.0);    // we init unused storage cells with NaN - why is it -NaN?

  // matrix A in band storage for gbsv - the banded format is N times 2*kl+ku+1 which is 
  // 9 times 2*3+2+1 = 9 times 8 in thsi case - because LAPACK uses in general, column-major 
  // indexing, our matrix looks transposed to the one in the comment above:
  double ab[ldab*N] = 
  {  _, _, _, _, _,11,21,31,41,     // 1st column of banded format
     _, _, _, _,12,22,32,42,52,     // 2nd column
     _, _, _,13,23,33,43,53,63,
     _, _, _,24,34,44,54,64,74,
     _, _, _,35,45,55,65,75,85,
     _, _, _,46,56,66,76,86,96,
     _, _, _,57,67,77,87,97, _,
     _, _, _,68,78,88,98, _, _,
     _, _, _,79,89,99, _, _, _ };  // 9th column

  // the matrix X in A * X = B:
  double X[ldb*nrhs] = {
    1,2,3,4,5,6,7,8,9,
    9,8,7,6,5,4,3,2,1,
    1,0,2,0,3,0,4,0,5,
    0,1,0,2,0,3,0,4,0
  };

  // right-hand-side B in A * X = B, after solving the system, it will comtain the solution, so it
  // should be equal to X up to roundoff errors:
  double B[ldb*nrhs];

  // let's do a vector-rhs first (instead of a matrix rhs):
  double x[ldb] = { 1,2,3,4,5,6,7,8,9 };
  double b[ldb] = { _,_,_,_,_,_,_,_,_ };  // the right hand side in A * x = b

  // to obtain the target values for b:
  //for(int i = 0; i < N; i++) {
  //  b[i] = 0;
  //  for(int j = 0; j < N; j++)
  //    b[i] += refMat[j][i] * x[j];
  //}
  //// b is 74,230,505,931,1489,2179,3001,3055,2930



  // for the banded storage used by gbmv, we need to store the matrix in the format:

  // ** ** 13 24 35 46 57 68 79   upper diagonal 2
  // ** 12 23 34 45 56 67 78 89   upper diagonal 1
  // 11 22 33 44 55 66 77 88 99   main diagonal
  // 21 32 43 54 65 76 87 98 **   lower diagonal 1
  // 31 42 53 64 75 86 97 ** **   lower diagonal 2
  // 41 52 63 74 85 96 ** ** **   lower diagonal 3

  static Int lda = kl+ku+1; // leading dimension of a >= kl+ku+1 
  int M = N;                // number of rows
  double alpha = 1.0;       // scaler in gbmv
  double beta  = 0.0;       // gbmv computes Y = alpha*A*x + beta*y

  // matrix A in band storage for gbmv:
  double a[lda*N] =
  {  _, _,11,21,31,41,     // 1st column of banded format
     _,12,22,32,42,52,     // 2nd column
    13,23,33,43,53,63,
    24,34,44,54,64,74,
    35,45,55,65,75,85,
    46,56,66,76,86,96,
    57,67,77,87,97, _,
    68,78,88,98, _, _,
    79,89,99, _, _, _ };   // 9th column

  //double one = 1, zero = 0;
  long int incX = 1;         // might be wrong
  long int incY = 1;         //
  //double yDummy;             // dummy - not referenced
  ftnlen trans_len = 0;      // ftnlen is typedef'd as "long" in f2c.h - i don't know, what this is
                             // used for, there's no documentation for that gbmv parameter -> check source 
                             // -> it's actually not used anywhere
  char trans = 'N';

  // gbmv needs pointers to "integer" which is defined as "long int" in f2c.h:
  long int N_ = N, M_ = M, lda_ = lda, kl_ = kl, ku_ = ku;
  gbmv(&trans, &M_, &N_, &kl_, &ku_, &alpha, a, &lda_, x, &incX, 
    &beta, b, &incY, trans_len);  

  // b should be 74,230,505,931,1489,2179,3001,3055,2930
  // but is actually 74,230,505,931,1489,2179,_,_,_  ...close, but not quite
  // the local variable leny is 6 in the routine - but should be 9 (length of y)
  // it gets leny from m which is the lda input - we should probably pass N_ for lda, leading 
  // dimension may not mean the array size but the logical dimension of the matrix a

  // aahh - that's wrong! that will work only if the rhs is a vector, but we want a matrix rhs here
  // ...so we actually need a gbmm - but such a routine doesn't seem to exist in BLAS - maybe i should 
  // write one myself as convenience function? or should we start with a vector rhs first? yes - 
  // that's probably better



  // can the upper version be passed to gbmv with some trickery, too? ...to avoid storing it twice?

  //int gbsv(long int *n, long int *kl, long int *ku, long int *nrhs, T *ab, long int *ldab, 
  //long int *ipiv, T *b, long int *ldb, long int *info);

  return r;
}

bool runUnitTests()
{
  bool r = true;

  r &= gbsvUnitTest();

  return r;
}


int main()
{  

  // move to blas1UnitTest or axpyUnitTest
  typedef std::vector<double> Vec;

  long int N = 10;
  //Vec a = rangeVector(N, 1, 1);
  double a = 2;
  Vec x = rangeVector(N, 2, 2);
  Vec y = rangeVector(N, 3, 3);

  //int daxpy_(integer *n, doublereal *da, doublereal *dx, 
  //  integer *incx, doublereal *dy, integer *incy)

  long int incX = 1, incY = 1;
  int result = daxpy_(&N, &a, &x[0], &incX, &y[0], &incY);
  // int result = daxpy_(&N, &a[0], &x[0], &incX, &y[0], &incY);

  result = axpy(&N, &a, &x[0], &incX, &y[0], &incY);


  // for level2, maybe use dgemv (dense general matrix-vector) and for level3 dgemm (dense general
  // matrix-matrix)

  // ok, this looks good
  // maybe make a function bool axpyUnitTest

  bool unitTestsPassed = runUnitTests();


  return 0;
}

