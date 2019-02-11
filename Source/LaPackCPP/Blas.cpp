#include "Blas.hpp"
namespace BlasCPP {

// Reference BLAS level1 routine (version 3.8.0)
template<class T>
int axpy(long int* n, T *da, T *dx, long int *incx, T *dy, long int *incy)
{
  // System generated locals
  long int i__1;

  // Local variables
  static long int i__, m, ix, iy, mp1;

  // Parameter adjustments
  --dy;
  --dx;

  // Function Body
  if(*n <= 0) {
    return 0;
  }
  if(*da == 0.) {
    return 0;
  }
  if(*incx == 1 && *incy == 1) {
  // code for both increments equal to 1
    m = *n % 4;
    if(m != 0) {
      i__1 = m;
      for(i__ = 1; i__ <= i__1; ++i__) {  // clean-up loop
        dy[i__] += *da * dx[i__];
      }
    }
    if(*n < 4) {
      return 0;
    }
    mp1 = m + 1;
    i__1 = *n;
    for(i__ = mp1; i__ <= i__1; i__ += 4) {
      dy[i__] += *da * dx[i__];
      dy[i__ + 1] += *da * dx[i__ + 1];
      dy[i__ + 2] += *da * dx[i__ + 2];
      dy[i__ + 3] += *da * dx[i__ + 3];
    }
  }
  else {
  // code for unequal increments or equal increments not equal to 1
    ix = 1;
    iy = 1;
    if(*incx < 0) {
      ix = (-(*n) + 1) * *incx + 1;
    }
    if(*incy < 0) {
      iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for(i__ = 1; i__ <= i__1; ++i__) {
      dy[iy] += *da * dx[ix];
      ix += *incx;
      iy += *incy;
    }
  }
  return 0;
}

//-------------------------------------------------------------------------------------------------

// Reference BLAS level2 routine (version 3.7.0)
template<class T>
int gbmv(char *trans, integer *m, integer *n, integer *kl, integer *ku, T *alpha, T *a, 
  integer *lda, T *x, integer *incx, T *beta, T *y, integer *incy, ftnlen trans_len)
{
  // System generated locals
  integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;

  // Local variables
  static integer i__, j, k, ix, iy, jx, jy, kx, ky, kup1, info;
  static T temp;
  static integer lenx, leny;
  extern logical lsame_(char *, char *, ftnlen, ftnlen);
  extern int xerbla_(char *, integer *, ftnlen);

  // Parameter adjustments
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --x;
  --y;

  // Function Body
  info = 0;
  if (! lsame_(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame_(trans, "T", (
    ftnlen)1, (ftnlen)1) && ! lsame_(trans, "C", (ftnlen)1, (ftnlen)1)
    ) {
    info = 1;
  } else if (*m < 0) {
    info = 2;
  } else if (*n < 0) {
    info = 3;
  } else if (*kl < 0) {
    info = 4;
  } else if (*ku < 0) {
    info = 5;
  } else if (*lda < *kl + *ku + 1) {
    info = 8;
  } else if (*incx == 0) {
    info = 10;
  } else if (*incy == 0) {
    info = 13;
  }
  if (info != 0) {
    xerbla_("DGBMV ", &info, (ftnlen)6);
    return 0;
  }

  // Quick return if possible.
  if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
    return 0;
  }

  // Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  // up the start points in  X  and  Y.
  if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
    lenx = *n;
    leny = *m;
  } else {
    lenx = *m;
    leny = *n;
  }
  if (*incx > 0) {
    kx = 1;
  } else {
    kx = 1 - (lenx - 1) * *incx;
  }
  if (*incy > 0) {
    ky = 1;
  } else {
    ky = 1 - (leny - 1) * *incy;
  }

  // Start the operations. In this version the elements of A are 
  // accessed sequentially with one pass through the band part of A. 
  // First form  y := beta*y.

  if (*beta != 1.) {
    if (*incy == 1) {
      if (*beta == 0.) {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[i__] = 0.;
          // L10: 
        }
      } else {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[i__] = *beta * y[i__];
          // L20: 
        }
      }
    } else {
      iy = ky;
      if (*beta == 0.) {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[iy] = 0.;
          iy += *incy;
          // L30: 
        }
      } else {
        i__1 = leny;
        for (i__ = 1; i__ <= i__1; ++i__) {
          y[iy] = *beta * y[iy];
          iy += *incy;
          // L40: 
        }
      }
    }
  }
  if (*alpha == 0.) {
    return 0;
  }
  kup1 = *ku + 1;
  if (lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {

    // Form  y := alpha*A*x + y.
    jx = kx;
    if (*incy == 1) {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = *alpha * x[jx];
        k = kup1 - j;
        // Computing MAX 
        i__2 = 1, i__3 = j - *ku;
        // Computing MIN 
        i__5 = *m, i__6 = j + *kl;
        i__4 = min(i__5,i__6);
        for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
          y[i__] += temp * a[k + i__ + j * a_dim1];
          // L50: 
        }
        jx += *incx;
        // L60: 
      }
    } else {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = *alpha * x[jx];
        iy = ky;
        k = kup1 - j;
        // Computing MAX 
        i__4 = 1, i__2 = j - *ku;
        // Computing MIN 
        i__5 = *m, i__6 = j + *kl;
        i__3 = min(i__5,i__6);
        for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
          y[iy] += temp * a[k + i__ + j * a_dim1];
          iy += *incy;
          // L70: 
        }
        jx += *incx;
        if (j > *ku) {
          ky += *incy;
        }
        // L80: 
      }
    }
  } else {

    // Form  y := alpha*A**T*x + y.
    jy = ky;
    if (*incx == 1) {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = 0.;
        k = kup1 - j;
        // Computing MAX 
        i__3 = 1, i__4 = j - *ku;
        // Computing MIN 
        i__5 = *m, i__6 = j + *kl;
        i__2 = min(i__5,i__6);
        for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
          temp += a[k + i__ + j * a_dim1] * x[i__];
          // L90:
        }
        y[jy] += *alpha * temp;
        jy += *incy;
        // L100: 
      }
    } else {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        temp = 0.;
        ix = kx;
        k = kup1 - j;
        // Computing MAX 
        i__2 = 1, i__3 = j - *ku;
        // Computing MIN 
        i__5 = *m, i__6 = j + *kl;
        i__4 = min(i__5,i__6);
        for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
          temp += a[k + i__ + j * a_dim1] * x[ix];
          ix += *incx;
          // L110:
        }
        y[jy] += *alpha * temp;
        jy += *incy;
        if (j > *ku) {
          kx += *incx;
        }
        // L120: 
      }
    }
  }

  return 0;

} // gbmv










// explicit template instantiations - todo: move them into separate files, one file per
// datatype:
template int axpy(long int* n, double *da, double *dx, long int *incx, 
  double *dy, long int *incy);

template int gbmv(char *trans, integer *m, integer *n, integer *kl, integer *ku, double *alpha, 
  double *a, integer *lda, double *x, integer *incx, double *beta, double *y, integer *incy, 
  ftnlen trans_len);








} // namespace BlasCPP