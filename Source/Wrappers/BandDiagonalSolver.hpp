#pragma once


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
  void setSystemSizes(int matrixOrder, int numSubDiagonals, int numSuperDiagonals, 
    int numRightHandSides);
  // maybe don't pass the number of right-hand sides here

  /** Sets one of the values in one of the digaonals. The diagIndex indicats which diagonal is 
  meant where 0 is the main diagonal, -1 is the first subdiagonal, +1 the first superdiagonal, 
  -2 the second subdiagonal and so on. Within each diagonal, the elemIndex counts from the 
  top-left to the bottom right, where the main-diagonal has N elements, the first sub- and 
  superdiagonals N-1 elements and so on. */
  void setDiagonalElement(int diagIndex, int elemIndex, T value);

  //void setAlgorithm(


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  // to be called after solve
  // getReciprocalPivotGrowth
  // getConditionNumber


  //-----------------------------------------------------------------------------------------------
  /** \name Computation */

  void solve(int numRightHandSides, const T* rightHandSides, T* solutions);

protected:

  // system size parameters:
  long N;       // matrix size
  long kl;      // number of subdiagonals
  long ku;      // number of superdiagonals
  long nrhs;    // number of right hand sides

  // workspace buffers:
  std::vector<T> a;    // the matrix in band-storage
  std::vector<T> b;    // right hand sides
  std::vector<T> afb;

};