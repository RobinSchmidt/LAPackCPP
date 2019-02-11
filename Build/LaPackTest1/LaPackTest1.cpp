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



  // for the banded storage used by gbmv, we need to store the matrix in the format:

  // ** ** 13 24 35 46 57 68 79   upper diagonal 2
  // ** 12 23 34 45 56 67 78 89   upper diagonal 1
  // 11 22 33 44 55 66 77 88 99   main diagonal
  // 21 32 43 54 65 76 87 98 **   lower diagonal 1
  // 31 42 53 64 75 86 97 ** **   lower diagonal 2
  // 41 52 63 74 85 96 ** ** **   lower diagonal 3


  // for the banded storage format used by gbsv, we to store the matrix in the form:

  // ** ** ** ** ** ++ ++ ++ ++   fields with ++ are used internally by the algorithm, fields
  // ** ** ** ** ++ ++ ++ ++ ++   with ** are not accessed
  // ** ** ** ++ ++ ++ ++ ++ ++
  // ** ** 13 24 35 46 57 68 79   upper diagonal 2
  // ** 12 23 34 45 56 67 78 89   upper diagonal 1
  // 11 22 33 44 55 66 77 88 99   main diagonal
  // 21 32 43 54 65 76 87 98 **   lower diagonal 1
  // 31 42 53 64 75 86 97 ** **   lower diagonal 2
  // 41 52 63 74 85 96 ** ** **   lower diagonal 3

  // can the lower version be passed to gbmv with some trickery, too? ...to avoid storing it twice?



  typedef const long int Int;  // for convenience
  static Int N    = 9;         // size of the system: NxN
  static Int kl   = 2;         // number of lower diagonals
  static Int ku   = 3;         // number of upper diagonals
  static Int nrhs = 4;         // number of right hand sides
  static Int ldab = 2*kl+ku+1; // leading dimension of ab >= 2*kl+ku+1
  static Int ldb  = N;         // leading dimension of b >= max(1,N). ...is N correct?
  long int ipiv[N];            // pivot indices
  double ab[ldab*N];           // matrix A in band storage for gbsv
  double b[ldb*nrhs];          // right-hand-side and solution matrix
  int info = 666;              // 0 on exit, if successful

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




  // todo: 
  // -fill the matrix a for gbmv using lapack's banded storage scheme
  // -fill the matrix ab for gbsv using lapack's banded storage scheme
  //  ...or can we pass the matrix ab also to gbmv...with some offsets, maybe?
  // -fill the right-hand-sides vectors b

  // we need gbmv: general banded matrix-vector multiply (blas2)

  static Int lda = kl+ku+1; // leading dimension of a >= kl+ku+1
  double a[lda*N];           // matrix A in band storage for gbmv
  int M = N;           // number of rows
  double alpha = 1.0;  // scaler in gbmv
  //int gbmv('N', &M, &N, &kl, &ku, &alpha, a, &lda, 
  // T *x, integer *incx, T *beta, T *y, integer *incy, ftnlen trans_len)


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

