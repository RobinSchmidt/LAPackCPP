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

  /** Sets one of the values in one of the digaonals. The diagIndex indicats which diagonal is 
  meant where 0 is the main diagonal, -1 is the first subdiagonal, +1 the first superdiagonal, 
  -2 the second subdiagonal and so on. Within each diagonal, the elemIndex counts from the 
  top-left to the bottom right, where the main-diagonal has N elements, the first sub- and 
  superdiagonals N-1 elements and so on. */
  void setDiagonalElement(int diagIndex, int elemIndex, T value);

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

  
  //-----------------------------------------------------------------------------------------------
  /** \name Computation */

  /** After the matrix has been set up via a call to setSystemSize and a bunch of calls to 
  setDiagonalElement, a call to solve will actually solve the system for a given number of right 
  hand sides and produce an equal number of solution vectors. */
  void solve(int numRightHandSides, T* rightHandSides, T* solutions);
  // todo: try to make rightHandSides const - figure out, if it's allowed that rightHandSides may
  // point to the same array as solutions

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