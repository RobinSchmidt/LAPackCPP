#pragma once

// may or may not be needed:
//#include "blas_extended.h"
//#include "blas_extended_private.h"





/**
Purpose:
gbmv computes y = alpha * A * x + beta * y, where 
A is a m x n banded matrix
x is a n x 1 vector
y is a m x 1 vector
alpha and beta are scalars 


Arguments:
order  Order of AP; row or column major
trans  Transpose of AB; no trans, trans, or conjugate trans
m      Dimension of AB
n      Dimension of AB and the length of vector x
kl     Number of lower diagnols of AB
ku     Number of upper diagnols of AB
alpha
AB
lda    Leading dimension of AB lda >= ku + kl + 1
x
incx   The stride for vector x.
beta
y
incy   The stride for vector y.
prec   Specifies the internal precision to be used.
        = blas_prec_single: single precision.
        = blas_prec_double: double precision.
        = blas_prec_extra : anything at least 1.5 times as accurate
                            than double, and wider than 80-bits.
                            We use double-double in our implementation. */
void BLAS_dgbmv_x(enum blas_order_type order, enum blas_trans_type trans, int m, int n, int kl, 
  int ku, double alpha, const double *a, int lda, const double *x, int incx, double beta, 
  double *y, int incy, enum blas_prec_type prec);