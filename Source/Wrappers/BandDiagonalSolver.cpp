// setup:

template<class T>
void rsBandDiagonalSolver<T>::setSystemSize(int matrixSize, int numSubDiagonals,
  int numSuperDiagonals)
{
  N    = matrixSize;
  kl   = numSubDiagonals;
  ku   = numSuperDiagonals;
  lda  = kl+ku+1;
  ldab = 2*kl+ku+1;
  ldb  = N;          // redundant - maybe get rid
  allocateMatrix();
}

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
  nrhs = numRightHandSides;

  allocateBuffers();
  if(algo == Algorithm::gbsv)
  {

  }
  else if(algo == Algorithm::gbsvx)
  {

  }
  else if(algo == Algorithm::gbsvxx) {
    gbsvxx(&fact, &trans, &N, &kl, &ku, &nrhs, &A[0], &lda, &AF[0], &ldab, &ipiv[0], &equed, &R[0], 
      &C[0], &B[0], &ldb, &X[0], &ldb, &rcond, &rpvgrw, &berr[0], &n_err_bnds, &err_bnds_norm[0], 
      &err_bnds_comp[0], &nparams, &params[0], &work[0], &iwork[0], &info, 0, 0, 0); 
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
  // todo: maybe re-allocate only those buffers that are needed for the selected routine - gbsv 
  // needs less buffers than gbsvxx, for example
                                          // used in routines (verify, if this is complete):
  ipiv.resize(N);                         // all
  R.resize(N);                            // gbsvx, gbsvxx 
  C.resize(N);                            // gbsvx, gbsvxx 
  ferr.resize(nrhs);                      // gbsvx only (i think)
  berr.resize(nrhs);                      // gbsvx, gbsvxx
  err_bnds_norm.resize(n_err_bnds*nrhs);  // gbsvxx
  err_bnds_comp.resize(n_err_bnds*nrhs);  // gbsvxx
  work.resize(4*N);                       // gbsvxx: 4*N, gbsvx: 3*N
  iwork.resize(N);                        // gbsvx, gbsvxx
}
