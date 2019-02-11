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

//  -- Reference BLAS level1 routine (version 3.1) -- 
logical lsame(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len)
{
  // System generated locals
  logical ret_val;

  // Local variables 
  static integer inta, intb, zcode;

  // Test if the characters are equal
  ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
  if(ret_val) {
    return ret_val;
  }

  // Now test for equivalence if both characters are alphabetic.
  zcode = 'Z';

  // Use 'Z' rather than 'A' so that ASCII can be detected on Prime 
  // machines, on which ICHAR returns a value with bit 8 set. 
  // ICHAR('A') on Prime machines returns 193 which is the same as 
  // ICHAR('A') on an EBCDIC machine.

  inta = *(unsigned char *)ca;
  intb = *(unsigned char *)cb;

  if(zcode == 90 || zcode == 122) {

    // ASCII is assumed - ZCODE is the ASCII code of either lower or upper case 'Z'. 
    if(inta >= 97 && inta <= 122) {
      inta += -32;
    }
    if(intb >= 97 && intb <= 122) {
      intb += -32;
    }

  }
  else if(zcode == 233 || zcode == 169) {

   // EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or 
   // upper case 'Z'.
    if(inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta
      >= 162 && inta <= 169) {
      inta += 64;
    }
    if(intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb
      >= 162 && intb <= 169) {
      intb += 64;
    }

  }
  else if(zcode == 218 || zcode == 250) {

   // ASCII is assumed, on Prime machines - ZCODE is the ASCII code
   // plus 128 of either lower or upper case 'Z'. 

    if(inta >= 225 && inta <= 250) {
      inta += -32;
    }
    if(intb >= 225 && intb <= 250) {
      intb += -32;
    }
  }
  ret_val = inta == intb;

  // End of LSAME

  return ret_val;
}

//-------------------------------------------------------------------------------------------------

// todo: fix linker errors for s_wsfe, s_stop, len_trim__, do_fio - these seem to be functions from
// libF2C
// -- Reference BLAS level1 routine (version 3.7.0) -- 
int xerbla(char *srname, integer *info, ftnlen srname_len)
{
  // Code of the function has been commented out by Robin Schmidt - at the moment, xerbla is just
  // a dummy function - i should probably set a debug-breakpoint here and re-implement it 
  // completely - it deals with some weird I/O functions from the f2c library which seems to be 
  // useless clutter...

  /*
  // Table of constant values
  static integer c__1 = 1;

  // Format strings 
  static char fmt_9999[] = "(\002 ** On entry to \002,a,\002 parameter num"
    "ber \002,i2,\002 had \002,\002an illegal value\002)";

  // Builtin functions 
  integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
  int s_stop(char *, ftnlen);

  // Local variables
  extern integer len_trim__(char *, ftnlen);
  // Note by Robin Schmidt:
  // this was declared as an intrinsic function in the original xerbla.f file, but it stumped the 
  // f2c translator, giving an error about an unknown intrinsic function, so i changed the 
  // "INTRINSIC" keyword to "EXTERNAL" - that allowed the translation, but i'm not sure, if it 
  // gives the correct behavior 

  // Fortran I/O blocks
  static cilist io___1 = { 0, 6, 0, fmt_9999, 0 };

  s_wsfe(&io___1);
  do_fio(&c__1, srname, len_trim__(srname, srname_len));
  do_fio(&c__1, (char *)&(*info), (ftnlen)sizeof(integer));
  e_wsfe();

  s_stop("", (ftnlen)0);

  // End of XERBLA
  */

  return 0;
} // xerbla

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
  extern logical lsame(char *, char *, ftnlen, ftnlen);
  extern int xerbla(char *, integer *, ftnlen);

  // Parameter adjustments
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --x;
  --y;

  // Function Body
  info = 0;
  if (! lsame(trans, "N", (ftnlen)1, (ftnlen)1) && ! lsame(trans, "T", (
    ftnlen)1, (ftnlen)1) && ! lsame(trans, "C", (ftnlen)1, (ftnlen)1)
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
    xerbla("DGBMV ", &info, (ftnlen)6);
    return 0;
  }

  // Quick return if possible.
  if (*m == 0 || *n == 0 || *alpha == 0. && *beta == 1.) {
    return 0;
  }

  // Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  // up the start points in  X  and  Y.
  if (lsame(trans, "N", (ftnlen)1, (ftnlen)1)) {
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
  if (lsame(trans, "N", (ftnlen)1, (ftnlen)1)) {

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
// datatype - maybe the particular type for which we instantiate the templates can be #defined
// the i just need to write the instantiation code once and copy/rename the files and just change 
// the #define - easier to maintain

template int axpy(long int* n, double *da, double *dx, long int *incx, 
  double *dy, long int *incy);

template int gbmv(char *trans, integer *m, integer *n, integer *kl, integer *ku, double *alpha, 
  double *a, integer *lda, double *x, integer *incx, double *beta, double *y, integer *incy, 
  ftnlen trans_len);








} // namespace BlasCPP