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
  rsBandDiagonalSolver(int matrixSize, int numRightHandSides);

  void setMatrixSizeAndNumRightHandSides(int matrixSize, int numRightHandSides);

  //void setAlgorithm(


protected:

  long N;            // matrix size
  long nrhs;         // number of right hand sides

  std::vector<T> a;  // the matrix in band-storage
  std::vector<T> b;  // right hand sides

};