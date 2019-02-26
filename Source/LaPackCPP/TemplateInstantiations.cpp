#pragma once

namespace LaPackCPP {


//=================================================================================================

// Blas

// explicit template instantiations - todo: move them into separate files, one file per
// datatype - maybe the particular type for which we instantiate the templates can be #defined
// the i just need to write the instantiation code once and copy/rename the files and just change 
// the #define - easier to maintain

// Level 1:

template int axpy(long int* n, double *da, double *dx, long int *incx,
  double *dy, long int *incy);

template double asum(integer *n, double *dx, integer *incx);

template int copy(long *n, double *dx, long *incx, double *dy, long *incy);

template double dot(integer *n, double *dx, integer *incx, double *dy, integer *incy);


template long iamax(long *n, double *dx, long *incx);

template int scal(long *n, double *da, double *dx, long *incx);

template int swap(long *n, double *dx, long *incx, double *dy, long *incy);

// Level 2:

//template int ger(integer *m, integer *n, doublereal *alpha, doublereal *x, integer *incx, 
//  doublereal *y, integer *incy, doublereal *a, integer *lda);

template int ger(long *m, long *n, double *alpha, double *x, long *incx, double *y, long *incy,
  double *a, long *lda);

template int gbmv(char *trans, integer *m, integer *n, integer *kl, integer *ku, double *alpha,
  double *a, integer *lda, double *x, integer *incx, double *beta, double *y, integer *incy,
  ftnlen trans_len);

template int gemv(char *trans, integer *m, integer *n, double *alpha, double *a, integer *lda,
  double *x, integer *incx, double *beta, double *y, integer *incy, ftnlen trans_len);

template int tbsv(char *uplo, char *trans, char *diag, integer *n, integer *k, double *a,
  integer *lda, double *x, integer *incx, ftnlen uplo_len, ftnlen trans_len, ftnlen diag_len);


// Level 3:

template int gemm(char *transa, char *transb, integer *m, integer *n, integer *k, double *alpha,
  double *a, integer *lda, double *b, integer *ldb, double *beta, double *c__, integer *ldc,
  ftnlen transa_len, ftnlen transb_len);

template int trsm(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n,
  double *alpha, double *a, integer *lda, double *b, integer *ldb, ftnlen side_len,
  ftnlen uplo_len, ftnlen transa_len, ftnlen diag_len);

}