#pragma once

#include <vector>

/** A convenience wrapper class for LAPACK's routines for solving band-diagonal systems of linear 
equations. It wraps the routines _gbsv, _gbsvx and _gbsvxx where the underscore is to be understood
as a placeholdder for the datatype, i.e. s,d,c,z for single/double/single-complex/double-complex 
numbers. So far, it has been tested only for double real numbers...todo: check, if it actually 
works for all these 4 datatypes... */

template<class T>
class rsBandDiagonalSolver
{

public:

  /** Enumeration of the available algorithms. */
  enum class Algorithm
  {
    gbsv,    // LAPACK's simplest driver routine for band diagonal systems
    gbsvx,   // LAPACK's expert driver routine
    gbsvxx   // LAPACK's extended expert driver 
  };

  /** Constructor. You have to pass the size of the system, i.e. the number of rows and columns of
  the matrix. Allocates memory for the temporary workspace arrays. */
  //rsBandDiagonalSolver(int matrixSize, int numRightHandSides);

  //void setMatrixSizeAndNumRightHandSides(int matrixSize, int numRightHandSides);

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the various size parameters of the system to be solved. The matrixOrder is the number
  of rows and columns of the matrix, the others are self-explanatory */
  void setSystemSize(int matrixOrder, int numSubDiagonals, int numSuperDiagonals);
    //int numRightHandSides);
  // maybe don't pass the number of right-hand sides here

  /** Sets one of the values in one of the diagonals. The diagIndex indicates which diagonal is 
  meant where 0 is the main diagonal, -1 is the first subdiagonal, +1 the first superdiagonal, 
  -2 the second subdiagonal and so on. Within each diagonal, the elemIndex counts from the 
  top-left to the bottom right, where the main-diagonal has N elements, the first sub- and 
  superdiagonals N-1 elements and so on. */
  void setDiagonalElement(int diagIndex, int elemIndex, const T& value)
  {
    A[diagElemIndex(diagIndex, elemIndex)] = value;
  }

  /** not yet implemented
  Accesses matrix elements via regular row- and column indices i, j. Raises an error, if the 
  index pair is outside the band of nonzero values */
  void setElement(int i, int j, const T& value)
  {

  }
  // 
  // AB(KL+KU+1+i-j,j) = A(i,j)
  /*
  storage scheme for gbmv (used for gbsvx and gbsvxx):
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


  storage scheme for gbsv:
  AB:   array, dimension (LDAB,N). On entry, the matrix A in band storage, in rows KL+1 to 
  2*KL+KU+1; rows 1 to KL of the array need not be set. The j-th column of A is stored in the 
  j-th column of the array AB as follows: AB(KL+KU+1+i-j,j) = A(i,j) for 
  max(1,j-KU)<=i<=min(N,j+KL). On exit, details of the factorization: U is stored as an upper 
  triangular band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and the multipliers 
  used during the factorization are stored in rows KL+KU+2 to 2*KL+KU+1. See below for further
  details.
  */

  //void setAlgorithm(


  /** Select whether or not equilibration should be used (if necessarry). */
  void setEquilibration(bool shouldEquilibrate);

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns reciprocal condition number of the matrix. Call it after solve(). */
  T getReciprocalConditionNumber() { return rcond; }

  /** Returns reciprocal pivot growth factor. Call it after solve(). */
  T getReciprocalPivotGrowth() { return rpvgrw; }

  //T getMaxErrorBound()
  // getRowScalers, getColumnScalers,....

  // todo: error-bounds, etc.

  // getDiagonalElement, getElement, operator()(int i, int j)

  
  //-----------------------------------------------------------------------------------------------
  /** \name Computation */

  /** After the matrix has been set up via a call to setSystemSize and a bunch of calls to 
  setDiagonalElement, a call to solve will actually solve the system for a given number of right 
  hand sides and produce an equal number of solution vectors. */
  void solve(int numRightHandSides, T* rightHandSides, T* solutions);
  // todo: try to make rightHandSides const - figure out, if it's allowed that rightHandSides may
  // point to the same array as solutions

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Conversion for index of the diagonal (ranging from -kl to ku) and index of element within the 
  diagonal (ranging from 0 to N-abs(diagIndex)) to the flat array index for LAPACK's band storage
  format. Mainly for internal use in setDiagonalElement, but you can also use it yourself - but 
  only after setting up setSystemSize appropriately, of course. */
  inline int diagElemIndex(int diagIndex, int elemIndex)
  {
    int d = diagIndex;
    int e = elemIndex;
    int row = ku - d; 
    int col = e;
    if(d > 0) 
      col += d;
    return bandedRowColToIndex(row, col, kl, ku);

    //int i = (kl+ku+1)*col + row;
    //// maybe raise an error if(i < 0 || i >= A.size())
    //return i;
  }
  // maybe factor this out into a class BandDiagonalMatrixView or something - maybe as static
  // function that takes kl, ku as additional arguments
  // maybe rename to diagIndicesToFlatArrayIndex, diagElemToArrayIndex, maybe move implementation
  // out of class

  /** Converts from regular (dense) matrix indices for the row and column to the actual index in 
  the storage array. */
  inline int rowColToArrayIndex(int row, int col)
  {
    //row += ku-col;               // row index must be manipulated according to the column
    //int i = (kl+ku+1)*col + row; // convert band-matrix indices to flat array index
    //// maybe raise an error if(i < 0 || i >= A.size())
    //return i;

    return bandedRowColToIndex(row+ku-col, col, kl, ku);
  }
  // the computation of i from the row- and column indices of the banded storage format is the same
  // in rowColToArrayIndex and diagElemIndex - maybe factor out into function bandedRowColToIndex
  // and call this function denseRowColToIndex


  /** Converts from a row- and column index given in banded storage format to a flat array 
  index. */
  static inline int bandedRowColToIndex(int row, int col, int numSubDiags, int numSuperDiags)
  {
    return (numSubDiags+numSuperDiags+1)*col + row;
  }


protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Memory allocation */

  /** Allocates the storage memory for the matrix. Called from setSystemSize(). */
  void allocateMatrix();

  /** Allocates the temporary workspace buffers. Called from solve(). */
  void allocateBuffers();

  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  // system size parameters:
  long N;       // matrix size
  long kl;      // number of subdiagonals
  long ku;      // number of superdiagonals
  long lda;     // leading dimension of matrix A
  long ldab;    // ?

  // right-hand-side size parameters:
  long nrhs;    // number of right hand sides
  long ldb;     // = N, maybe get rid

  // workspace buffers:
  std::vector<T> A;        // coefficient matrix in band-storage
  //std::vector<T> B;        // matrix of right hand sides
  std::vector<T> AF;       // factored form of A
  std::vector<T> work;     // workspace, size 4*N
  std::vector<T> R, C;     // row and column scale factors for equilibration, size N
  std::vector<long> iwork; // integer workspace, size N
  std::vector<long> ipiv;  // pivot indices, size N

  // options:
  Algorithm algo = Algorithm::gbsvxx; // extended expert driver - most precise
  char trans = 'N';        // input matrix is not transposed
  char fact  = 'E';        // input matrix is not factored and shall be equilibrated
  long nparams = 0;        // number of additional parameters for gbsvxx (not yet used)
  T params[1];             // dummy - not referenced, if nparams == 0

  // computed outputs:
  long info;                     // 0 on exit, if successful
  char equed;                    // returns the form of equlibration that was done
  T rcond;                       // reciprocal condition number
  T rpvgrw;                      // reciprocal pivot growth
  std::vector<T> ferr;           // componentwise relative forward error, size nrhs (verify)
  std::vector<T> berr;           // componentwise relative backward error, size nrhs (verify)
  long n_err_bnds = 3;           // number of error bounds
  std::vector<T> err_bnds_norm;  // various error bounds (up to 3 for each rhs), size 3*nrhs
  std::vector<T> err_bnds_comp;  // size 3*nrhs
};