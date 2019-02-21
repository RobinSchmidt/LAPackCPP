#include "BandDiagonalSolver.hpp"

// setup:

template<class T>
void rsBandDiagonalSolver<T>::setSystemSize(int matrixOrder, int numSubDiagonals,
  int numSuperDiagonals)
{
  N    = matrixOrder;
  kl   = numSubDiagonals;
  ku   = numSuperDiagonals;
  lda  = kl+ku+1;
  ldab = 2*kl+ku+1;
  allocateMatrix();
}

/*
template<class T>
void rsBandDiagonalSolver<T>::setDiagonalElement(int diagIndex, int elemIndex, T value)
{
  //A[diagElemIndex(diagIndex, elemIndex)] = value;

  int d = diagIndex;
  int e = elemIndex;
  int row = ku - d; 
  int col = e;
  if(d > 0) 
    col += d;
  int i = (kl+ku+1)*col + row;
  if(i < 0 || i >= A.size())
  {
    // debug-break/assert
    return;
  }
  A[i] = value;
}
*/

template<class T>
void rsBandDiagonalSolver<T>::setEquilibration(bool shouldEquilibrate)
{
  if(shouldEquilibrate)
    fact = 'E';
  else
    fact = 'N';
}

// inquiry:


// compute:

template<class T>
void rsBandDiagonalSolver<T>::solve(int numRightHandSides, T* B, T* X)
{
  // todo: set up nrhs, etc.
   
  nrhs = numRightHandSides;
  ldb  = N;  // maybe get rid
  allocateBuffers();



  if(algo == Algorithm::gbsv)
  {

  }
  else if(algo == Algorithm::gbsvx)
  {

  }
  else if(algo == Algorithm::gbsvxx)
  {
    gbsvxx(
      &fact, 
      &trans, 
      &N, 
      &kl, 
      &ku, 
      &nrhs, 
      &A[0], 
      &lda, 
      &AF[0], 
      &ldab,
      &ipiv[0], 
      &equed, 
      &R[0], 
      &C[0], 
      &B[0], 
      &ldb, 
      &X[0], 
      &ldb,          // check, if ldx == ldb? - should be
      &rcond, 
      &rpvgrw, 
      &berr[0], 
      &n_err_bnds, 
      &err_bnds_norm[0], 
      &err_bnds_comp[0], 
      &nparams,
      &params[0], 
      &work[0], 
      &iwork[0], 
      &info, 
      0, 0, 0); 
  }
  else
  {
    // raise error
  }
}

// memory allocation:

template<class T>
void rsBandDiagonalSolver<T>::allocateMatrix()
{
  A.resize(lda*N);
  AF.resize(ldab*N);
}

template<class T>
void rsBandDiagonalSolver<T>::allocateBuffers()
{

}
