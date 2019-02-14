
#include "Blas.hpp"
#include "LaPack.hpp"
#include <cmath>       // maybe move elsewhere

using namespace BlasCPP;

namespace LaPackCPP {

// some fiddling to make it compile and link - those functions are defined outside the LaPackCPP
// namespace, so in order to not have to enter a :: each time a function is called, we use this
// wrappers/delegators here inside the namspace - todo: clean this up! try to get rid and if
// impossible, at least inline them and maybe move to some other file

double log(doublereal x)  
{
  return ::log((double)x); 
}
integer i_len(char *s, ftnlen n)
{
  return ::i_len(s, n);
}
integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
{
  return ::s_cmp(a0, b0, la, lb);
}
int s_copy(register char *a, register char *b, ftnlen la, ftnlen lb)
{
  ::s_copy(a, b, la, lb);
  return 0;
}
integer i_nint(f2c_real *x)
{
  return ::i_nint(x);
}

//=================================================================================================

// translated from dgbsv, LAPACK driver routine (version 3.7.0) -- 
template<class T>
int gbsv(long int *n, long int *kl, long int *ku, long int *nrhs, T *ab, long int *ldab,
  long int *ipiv, T *b, long int *ldb, long int *info)
{
  // System generated locals
  long int ab_dim1, ab_offset, b_dim1, b_offset, i__1;

  // Local variables
  // Subroutine
  extern  int
    gbtrf(long int *, long int *, long int *,
      long int *, double *, long int *, long int *, long int *),
    xerbla_(char *, long int *, ftnlen),
    gbtrs(char *, long int *,
      long int *, long int *, long int *, double *, long int *, long int *, double *,
      long int *, long int *, ftnlen);

  // Parameter adjustments
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --ipiv;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;

  // Function Body
  *info = 0;
  if(*n < 0) {
    *info = -1;
  }
  else if(*kl < 0) {
    *info = -2;
  }
  else if(*ku < 0) {
    *info = -3;
  }
  else if(*nrhs < 0) {
    *info = -4;
  }
  else if(*ldab < (*kl << 1) + *ku + 1) {
    *info = -6;
  }
  else if(*ldb < max(*n, 1)) {
    *info = -9;
  }
  if(*info != 0) {
    i__1 = -(*info);
    xerbla_("DGBSV ", &i__1, (ftnlen)6);
    return 0;
  }

  // Compute the LU factorization of the band matrix A:
  gbtrf(n, n, kl, ku, &ab[ab_offset], ldab, &ipiv[1], info);
  if(*info == 0) {
    // Solve the system A*X = B, overwriting B with X:
    gbtrs("No transpose", n, kl, ku, nrhs, &ab[ab_offset], ldab, &ipiv[1], &b[b_offset], ldb,
      info, (ftnlen)12);
  }
  return 0;
}

//-------------------------------------------------------------------------------------------------

// translated from dgbtf2, LAPACK computational routine (version 3.7.0)
template<class T>
int gbtf2(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, integer *ipiv, 
  integer *info)
{
  /* Table of constant values */
  static integer c__1 = 1;
  static doublereal c_b9 = -1.;

  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4;
  T d__1;

  /* Local variables */
  static integer i__, j, km, jp, ju, kv;
  //extern /* Subroutine */ int dger_(integer *, integer *, T *, 
  //  T *, integer *, T *, integer *, T *, 
  //  integer *), dscal_(integer *, T *, T *, integer 
  //    *), dswap_(integer *, T *, integer *, T *, 
  //      integer *);
  //extern integer idamax_(integer *, T *, integer *);
  //extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);

  /*     KV is the number of superdiagonals in the factor U, allowing for */
  /*     fill-in. */

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --ipiv;

  /* Function Body */
  kv = *ku + *kl;

  /*     Test the input parameters. */

  *info = 0;
  if (*m < 0) {
    *info = -1;
  } else if (*n < 0) {
    *info = -2;
  } else if (*kl < 0) {
    *info = -3;
  } else if (*ku < 0) {
    *info = -4;
  } else if (*ldab < *kl + kv + 1) {
    *info = -6;
  }
  if (*info != 0) {
    i__1 = -(*info);
    xerbla("DGBTF2", &i__1, (ftnlen)6);
    return 0;
  }

  /*     Quick return if possible */

  if (*m == 0 || *n == 0) {
    return 0;
  }

  /*     Gaussian elimination with partial pivoting */

  /*     Set fill-in elements in columns KU+2 to KV to zero. */

  i__1 = min(kv,*n);
  for (j = *ku + 2; j <= i__1; ++j) {
    i__2 = *kl;
    for (i__ = kv - j + 2; i__ <= i__2; ++i__) {
      ab[i__ + j * ab_dim1] = 0.;
      /* L10: */
    }
    /* L20: */
  }

  /*     JU is the index of the last column affected by the current stage */
  /*     of the factorization. */

  ju = 1;

  i__1 = min(*m,*n);
  for (j = 1; j <= i__1; ++j) {

    /*        Set fill-in elements in column J+KV to zero. */

    if (j + kv <= *n) {
      i__2 = *kl;
      for (i__ = 1; i__ <= i__2; ++i__) {
        ab[i__ + (j + kv) * ab_dim1] = 0.;
        /* L30: */
      }
    }

    /*        Find pivot and test for singularity. KM is the number of */
    /*        subdiagonal elements in the current column. */

    /* Computing MIN */
    i__2 = *kl, i__3 = *m - j;
    km = min(i__2,i__3);
    i__2 = km + 1;
    jp = iamax(&i__2, &ab[kv + 1 + j * ab_dim1], &c__1);
    ipiv[j] = jp + j - 1;
    if (ab[kv + jp + j * ab_dim1] != 0.) {
      /* Computing MAX */
      /* Computing MIN */
      i__4 = j + *ku + jp - 1;
      i__2 = ju, i__3 = min(i__4,*n);
      ju = max(i__2,i__3);

      /*           Apply interchange to columns J to JU. */

      if (jp != 1) {
        i__2 = ju - j + 1;
        i__3 = *ldab - 1;
        i__4 = *ldab - 1;
        swap(&i__2, &ab[kv + jp + j * ab_dim1], &i__3, &ab[kv + 1 + 
          j * ab_dim1], &i__4);
      }

      if (km > 0) {

        /*              Compute multipliers. */

        d__1 = 1. / ab[kv + 1 + j * ab_dim1];
        scal(&km, &d__1, &ab[kv + 2 + j * ab_dim1], &c__1);

        /*              Update trailing submatrix within the band. */

        if (ju > j) {
          i__2 = ju - j;
          i__3 = *ldab - 1;
          i__4 = *ldab - 1;
          ger(&km, &i__2, &c_b9, &ab[kv + 2 + j * ab_dim1], &c__1,
            &ab[kv + (j + 1) * ab_dim1], &i__3, &ab[kv + 1 + 
            (j + 1) * ab_dim1], &i__4);
        }
      }
    } else {

      /*           If pivot is zero, set INFO to the index of the pivot */
      /*           unless a zero pivot has already been found. */

      if (*info == 0) {
        *info = j;
      }
    }
    /* L40: */
  }
  return 0;

  /*     End of DGBTF2 */

} /* gbtf2 */

//-------------------------------------------------------------------------------------------------

// translated from dgbtrf,LAPACK computational routine (version 3.7.0)
template<class T>
int gbtrf(integer *m, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, integer *ipiv,
  integer *info)
{
  // Table of constant values
  static integer c__1 = 1;
  static integer c__65 = 65;
  static T c_b18 = -1.;
  static T c_b31 = 1.;

  // System generated locals 
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
  T d__1;

  // Local variables 
  static integer i__, j, i2, i3, j2, j3, k2, jb, nb, ii, jj, jm, ip, jp, km,
    ju, kv, nw;
  // we need to comment these declarations - otherwise, the linker tries to find those functions
  // in the LaPackCPP namespace (but they belong to the BlasCPP namespace)
  //extern int ger(integer *, integer *, T *,
  //  T *, integer *, T *, integer *, T *,
  //  integer *);
  static T temp;
  //extern int scal(integer *, T *, T *,
  //  integer *), gemm(char *, char *, integer *, integer *, integer *
  //    , T *, T *, integer *, T *, integer *,
  //    T *, T *, integer *, ftnlen, ftnlen), dcopy_(
  //      integer *, T *, integer *, T *, integer *),
  //  swap(integer *, T *, integer *, T *, integer *
  //    );
  static T work13[4160], work31[4160];
  //extern int trsm(char *, char *, char *, char *,
  //  integer *, integer *, T *, T *, integer *,
  //  T *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dgbtf2_(
  //    integer *, integer *, integer *, integer *, T *, integer
  //    *, integer *, integer *);
  //extern integer iamax(integer *, T *, integer *);
  //extern int xerbla_(char *, integer *, ftnlen);
  //extern integer ilaenv(integer *, char *, char *, integer *, integer *,
  //  integer *, integer *, ftnlen, ftnlen);
  //extern int laswp(integer *, T *, integer *,
  //  integer *, integer *, integer *, integer *);


  // KV is the number of superdiagonals in the factor U, allowing for fill-in

  // Parameter adjustments
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --ipiv;

  // Function Body
  kv = *ku + *kl;

  // Test the input parameters.
  *info = 0;
  if(*m < 0) {
    *info = -1;
  }
  else if(*n < 0) {
    *info = -2;
  }
  else if(*kl < 0) {
    *info = -3;
  }
  else if(*ku < 0) {
    *info = -4;
  }
  else if(*ldab < *kl + kv + 1) {
    *info = -6;
  }
  if(*info != 0) {
    i__1 = -(*info);
    xerbla("DGBTRF", &i__1, (ftnlen)6);
    return 0;
  }

  // Quick return if possible
  if(*m == 0 || *n == 0) {
    return 0;
  }

  // Determine the block size for this environment
  nb = ilaenv(&c__1, "DGBTRF", " ", m, n, kl, ku, (ftnlen)6, (ftnlen)1);
  // The block size must not exceed the limit set by the size of the
  // local arrays WORK13 and WORK31.

  nb = min(nb, 64);

  if(nb <= 1 || nb > *kl) {
    // Use unblocked code
    gbtf2(m, n, kl, ku, &ab[ab_offset], ldab, &ipiv[1], info);
  }
  else {

  // Use blocked code
  // Zero the superdiagonal elements of the work array WORK13
    i__1 = nb;
    for(j = 1; j <= i__1; ++j) {
      i__2 = j - 1;
      for(i__ = 1; i__ <= i__2; ++i__) {
        work13[i__ + j * 65 - 66] = 0.;
        // L10: 
      }
      // L20: 
    }

    // Zero the subdiagonal elements of the work array WORK31
    i__1 = nb;
    for(j = 1; j <= i__1; ++j) {
      i__2 = nb;
      for(i__ = j + 1; i__ <= i__2; ++i__) {
        work31[i__ + j * 65 - 66] = 0.;
        // L30: 
      }
      // L40: 
    }

    // Gaussian elimination with partial pivoting
    // Set fill-in elements in columns KU+2 to KV to zero
    i__1 = min(kv, *n);
    for(j = *ku + 2; j <= i__1; ++j) {
      i__2 = *kl;
      for(i__ = kv - j + 2; i__ <= i__2; ++i__) {
        ab[i__ + j * ab_dim1] = 0.;
        // L50: 
      }
      // L60: 
    }

    // JU is the index of the last column affected by the current 
    // stage of the factorization 
    ju = 1;

    i__1 = min(*m, *n);
    i__2 = nb;
    for(j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2) {
      // Computing MIN 
      i__3 = nb, i__4 = min(*m, *n) - j + 1;
      jb = min(i__3, i__4);

      // The active part of the matrix is partitioned 
      // A11   A12   A13 
      // A21   A22   A23 
      // A31   A32   A33 
      // Here A11, A21 and A31 denote the current block of JB columns 
      // which is about to be factorized. The number of rows in the 
      // partitioning are JB, I2, I3 respectively, and the numbers 
      // of columns are JB, J2, J3. The superdiagonal elements of A13
      // and the subdiagonal elements of A31 lie outside the band.

      // Computing MIN 
      i__3 = *kl - jb, i__4 = *m - j - jb + 1;
      i2 = min(i__3, i__4);
      // Computing MIN
      i__3 = jb, i__4 = *m - j - *kl + 1;
      i3 = min(i__3, i__4);

      // J2 and J3 are computed after JU has been updated. 

      // Factorize the current block of JB columns 
      i__3 = j + jb - 1;
      for(jj = j; jj <= i__3; ++jj) {

        // Set fill-in elements in column JJ+KV to zero

        if(jj + kv <= *n) {
          i__4 = *kl;
          for(i__ = 1; i__ <= i__4; ++i__) {
            ab[i__ + (jj + kv) * ab_dim1] = 0.;
            // L70: 
          }
        }

        // Find pivot and test for singularity. KM is the number of 
        // subdiagonal elements in the current column. */

        // Computing MIN 
        i__4 = *kl, i__5 = *m - jj;
        km = min(i__4, i__5);
        i__4 = km + 1;
        jp = iamax(&i__4, &ab[kv + 1 + jj * ab_dim1], &c__1);
        ipiv[jj] = jp + jj - j;
        if(ab[kv + jp + jj * ab_dim1] != 0.) {
          // Computing MAX 
          // Computing MIN
          i__6 = jj + *ku + jp - 1;
          i__4 = ju, i__5 = min(i__6, *n);
          ju = max(i__4, i__5);
          if(jp != 1) {

            // Apply interchange to columns J to J+JB-1
            if(jp + jj - 1 < j + *kl) {

              i__4 = *ldab - 1;
              i__5 = *ldab - 1;
              swap(&jb, &ab[kv + 1 + jj - j + j * ab_dim1], &
                i__4, &ab[kv + jp + jj - j + j * ab_dim1],
                &i__5);
            }
            else {

            // The interchange affects columns J to JJ-1 of A31 
            // which are stored in the work array WORK31 
              i__4 = jj - j;
              i__5 = *ldab - 1;
              swap(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1],
                &i__5, &work31[jp + jj - j - *kl - 1], &
                c__65);
              i__4 = j + jb - jj;
              i__5 = *ldab - 1;
              i__6 = *ldab - 1;
              swap(&i__4, &ab[kv + 1 + jj * ab_dim1], &i__5, &
                ab[kv + jp + jj * ab_dim1], &i__6);
            }
          }

          // Compute multipliers 

          d__1 = 1. / ab[kv + 1 + jj * ab_dim1];
          scal(&km, &d__1, &ab[kv + 2 + jj * ab_dim1], &c__1);

          // Update trailing submatrix within the band and within 
          // the current block. JM is the index of the last column 
          // which needs to be updated. 

          // Computing MIN
          i__4 = ju, i__5 = j + jb - 1;
          jm = min(i__4, i__5);
          if(jm > jj) {
            i__4 = jm - jj;
            i__5 = *ldab - 1;
            i__6 = *ldab - 1;
            ger(&km, &i__4, &c_b18, &ab[kv + 2 + jj * ab_dim1],
              &c__1, &ab[kv + (jj + 1) * ab_dim1], &i__5, &
              ab[kv + 1 + (jj + 1) * ab_dim1], &i__6);
          }
        }
        else {

       // If pivot is zero, set INFO to the index of the pivot 
       // unless a zero pivot has already been found.

          if(*info == 0) {
            *info = jj;
          }
        }

        // Copy current column of A31 into the work array WORK31 

        // Computing MIN 
        i__4 = jj - j + 1;
        nw = min(i__4, i3);
        if(nw > 0) {
          copy(&nw, &ab[kv + *kl + 1 - jj + j + jj * ab_dim1], &
            c__1, &work31[(jj - j + 1) * 65 - 65], &c__1);
        }
        // L80: 
      }
      if(j + jb <= *n) {

        // Apply the row interchanges to the other blocks.

        // Computing MIN 
        i__3 = ju - j + 1;
        j2 = min(i__3, kv) - jb;
        // Computing MAX
        i__3 = 0, i__4 = ju - j - kv + 1;
        j3 = max(i__3, i__4);

        // Use DLASWP to apply the row interchanges to A12, A22, and 
        // A32.
        i__3 = *ldab - 1;
        laswp(&j2, &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__3, &
          c__1, &jb, &ipiv[j], &c__1);

        // Adjust the pivot indices.
        i__3 = j + jb - 1;
        for(i__ = j; i__ <= i__3; ++i__) {
          ipiv[i__] = ipiv[i__] + j - 1;
          // L90:
        }

        // Apply the row interchanges to A13, A23, and A33 
        // columnwise. 
        k2 = j - 1 + jb + j2;
        i__3 = j3;
        for(i__ = 1; i__ <= i__3; ++i__) {
          jj = k2 + i__;
          i__4 = j + jb - 1;
          for(ii = j + i__ - 1; ii <= i__4; ++ii) {
            ip = ipiv[ii];
            if(ip != ii) {
              temp = ab[kv + 1 + ii - jj + jj * ab_dim1];
              ab[kv + 1 + ii - jj + jj * ab_dim1] = ab[kv + 1 +
                ip - jj + jj * ab_dim1];
              ab[kv + 1 + ip - jj + jj * ab_dim1] = temp;
            }
            // L100: 
          }
          // L110: 
        }

        // Update the relevant part of the trailing submatrix 

        if(j2 > 0) {

          // Update A12 
          i__3 = *ldab - 1;
          i__4 = *ldab - 1;
          trsm("Left", "Lower", "No transpose", "Unit", &jb, &j2,
            &c_b31, &ab[kv + 1 + j * ab_dim1], &i__3, &ab[kv
            + 1 - jb + (j + jb) * ab_dim1], &i__4, (ftnlen)4,
            (ftnlen)5, (ftnlen)12, (ftnlen)4);

          if(i2 > 0) {

            // Update A22
            i__3 = *ldab - 1;
            i__4 = *ldab - 1;
            i__5 = *ldab - 1;
            gemm("No transpose", "No transpose", &i2, &j2, &jb,
              &c_b18, &ab[kv + 1 + jb + j * ab_dim1], &i__3,
              &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__4,
              &c_b31, &ab[kv + 1 + (j + jb) * ab_dim1], &
              i__5, (ftnlen)12, (ftnlen)12);
          }

          if(i3 > 0) {

            // Update A32
            i__3 = *ldab - 1;
            i__4 = *ldab - 1;
            gemm("No transpose", "No transpose", &i3, &j2, &jb,
              &c_b18, work31, &c__65, &ab[kv + 1 - jb + (j
                + jb) * ab_dim1], &i__3, &c_b31, &ab[kv + *kl
              + 1 - jb + (j + jb) * ab_dim1], &i__4, (
                ftnlen)12, (ftnlen)12);
          }
        }

        if(j3 > 0) {

          // Copy the lower triangle of A13 into the work array
          // WORK13
          i__3 = j3;
          for(jj = 1; jj <= i__3; ++jj) {
            i__4 = jb;
            for(ii = jj; ii <= i__4; ++ii) {
              work13[ii + jj * 65 - 66] = ab[ii - jj + 1 + (jj
                + j + kv - 1) * ab_dim1];
              // L120: 
            }
            // L130: 
          }

          // Update A13 in the work array 
          i__3 = *ldab - 1;
          trsm("Left", "Lower", "No transpose", "Unit", &jb, &j3,
            &c_b31, &ab[kv + 1 + j * ab_dim1], &i__3, work13,
            &c__65, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)
            4);

          if(i2 > 0) {
            // Update A23 
            i__3 = *ldab - 1;
            i__4 = *ldab - 1;
            gemm("No transpose", "No transpose", &i2, &j3, &jb,
              &c_b18, &ab[kv + 1 + jb + j * ab_dim1], &i__3,
              work13, &c__65, &c_b31, &ab[jb + 1 + (j + kv)
              * ab_dim1], &i__4, (ftnlen)12, (ftnlen)12);
          }

          if(i3 > 0) {
            // Update A33 
            i__3 = *ldab - 1;
            gemm("No transpose", "No transpose", &i3, &j3, &jb,
              &c_b18, work31, &c__65, work13, &c__65, &
              c_b31, &ab[*kl + 1 + (j + kv) * ab_dim1], &
              i__3, (ftnlen)12, (ftnlen)12);
          }

          // Copy the lower triangle of A13 back into place
          i__3 = j3;
          for(jj = 1; jj <= i__3; ++jj) {
            i__4 = jb;
            for(ii = jj; ii <= i__4; ++ii) {
              ab[ii - jj + 1 + (jj + j + kv - 1) * ab_dim1] =
                work13[ii + jj * 65 - 66];
              // L140: 
            }
            // L150: 
          }
        }
      }
      else {

      // Adjust the pivot indices.

        i__3 = j + jb - 1;
        for(i__ = j; i__ <= i__3; ++i__) {
          ipiv[i__] = ipiv[i__] + j - 1;
          // L160:
        }
      }

      // Partially undo the interchanges in the current block to 
      // restore the upper triangular form of A31 and copy the upper 
      // triangle of A31 back into place 
      i__3 = j;
      for(jj = j + jb - 1; jj >= i__3; --jj) {
        jp = ipiv[jj] - jj + 1;
        if(jp != 1) {

          // Apply interchange to columns J to JJ-1

          if(jp + jj - 1 < j + *kl) {
            // The interchange does not affect A31
            i__4 = jj - j;
            i__5 = *ldab - 1;
            i__6 = *ldab - 1;
            swap(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
              i__5, &ab[kv + jp + jj - j + j * ab_dim1], &
              i__6);
          }
          else {
            // The interchange does affect A31
            i__4 = jj - j;
            i__5 = *ldab - 1;
            swap(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
              i__5, &work31[jp + jj - j - *kl - 1], &c__65);
          }
        }

        // Copy the current column of A31 back into place
        // Computing MIN 
        i__4 = i3, i__5 = jj - j + 1;
        nw = min(i__4, i__5);
        if(nw > 0) {
          copy(&nw, &work31[(jj - j + 1) * 65 - 65], &c__1, &ab[
            kv + *kl + 1 - jj + j + jj * ab_dim1], &c__1);
        }
        // L170:
      }
      // L180:
    }
  }

  return 0;
} // gbtrf

//-------------------------------------------------------------------------------------------------

// translated from dgbtrs, LAPACK computational routine (version 3.7.0)
template<class T>
int gbtrs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab,
  integer *ldab, integer *ipiv, T *b, integer *ldb, integer *info, ftnlen trans_len)
{
  // Table of constant values
  static T c_b7 = -1.;
  static integer c__1 = 1;
  static T c_b23 = 1.;

  // System generated locals
  integer ab_dim1, ab_offset, b_dim1, b_offset, i__1, i__2, i__3;

  // Local variables 
  static integer i__, j, l, kd, lm;
  //extern int dger_(integer *, integer *, doublereal *,
  //  doublereal *, integer *, doublereal *, integer *, doublereal *,
  //  integer *);
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  //extern int dgemv_(char *, integer *, integer *,
  //  doublereal *, doublereal *, integer *, doublereal *, integer *,
  //  doublereal *, doublereal *, integer *, ftnlen), dswap_(integer *,
  //    doublereal *, integer *, doublereal *, integer *), dtbsv_(char *,
  //      char *, char *, integer *, integer *, doublereal *, integer *,
  //      doublereal *, integer *, ftnlen, ftnlen, ftnlen);
  static logical lnoti;
  //extern int xerbla(char *, integer *, ftnlen);
  static logical notran;

  // Parameter adjustments
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --ipiv;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;

  // Function Body
  *info = 0;
  notran = lsame(trans, "N", (ftnlen)1, (ftnlen)1);
  if(!notran && !lsame(trans, "T", (ftnlen)1, (ftnlen)1) && !lsame(
    trans, "C", (ftnlen)1, (ftnlen)1)) {
    *info = -1;
  }
  else if(*n < 0) {
    *info = -2;
  }
  else if(*kl < 0) {
    *info = -3;
  }
  else if(*ku < 0) {
    *info = -4;
  }
  else if(*nrhs < 0) {
    *info = -5;
  }
  else if(*ldab < (*kl << 1) + *ku + 1) {
    *info = -7;
  }
  else if(*ldb < max(1, *n)) {
    *info = -10;
  }
  if(*info != 0) {
    i__1 = -(*info);
    xerbla("DGBTRS", &i__1, (ftnlen)6);
    return 0;
  }

  // Quick return if possible
  if(*n == 0 || *nrhs == 0) {
    return 0;
  }

  kd = *ku + *kl + 1;
  lnoti = *kl > 0;

  if(notran) {

    // Solve  A*X = B. 
    // Solve L*X = B, overwriting B with X. 
    // L is represented as a product of permutations and unit lower 
    // triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1), 
    // where each transformation L(i) is a rank-one modification of 
    // the identity matrix. 
    if(lnoti) {
      i__1 = *n - 1;
      for(j = 1; j <= i__1; ++j) {
        // Computing MIN 
        i__2 = *kl, i__3 = *n - j;
        lm = min(i__2, i__3);
        l = ipiv[j];
        if(l != j) {
          swap(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
        }
        ger(&lm, nrhs, &c_b7, &ab[kd + 1 + j * ab_dim1], &c__1, &b[
          j + b_dim1], ldb, &b[j + 1 + b_dim1], ldb);
        // L10:
      }
    }

    i__1 = *nrhs;
    for(i__ = 1; i__ <= i__1; ++i__) {

      // Solve U*X = B, overwriting B with X.
      i__2 = *kl + *ku;
      tbsv("Upper", "No transpose", "Non-unit", n, &i__2, &ab[
        ab_offset], ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5,
          (ftnlen)12, (ftnlen)8);
      // L20: 
    }

  }
  else {

 // Solve A**T*X = B.
    i__1 = *nrhs;
    for(i__ = 1; i__ <= i__1; ++i__) {

      // Solve U**T*X = B, overwriting B with X.
      i__2 = *kl + *ku;
      tbsv("Upper", "Transpose", "Non-unit", n, &i__2, &ab[ab_offset],
        ldab, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)5, (ftnlen)9,
        (ftnlen)8);
      // L30:
    }

    // Solve L**T * X = B, overwriting B with X.
    if(lnoti) {
      for(j = *n - 1; j >= 1; --j) {
        // Computing MIN 
        i__1 = *kl, i__2 = *n - j;
        lm = min(i__1, i__2);
        gemv("Transpose", &lm, nrhs, &c_b7, &b[j + 1 + b_dim1], ldb,
          &ab[kd + 1 + j * ab_dim1], &c__1, &c_b23, &b[j +
          b_dim1], ldb, (ftnlen)9);
        l = ipiv[j];
        if(l != j) {
          swap(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
        }
        // L40:
      }
    }
  }
  return 0;

  // End of DGBTRS 

} // gbtrs

//=================================================================================================
// Auxiliary routines:


//-------------------------------------------------------------------------------------------------

// LAPACK auxiliary routine (version 3.7.0) 
integer ieeeck(integer *ispec, f2c_real *zero, f2c_real *one)
{
  /* System generated locals */
  integer ret_val;

  /* Local variables */
  static f2c_real nan1, nan2, nan3, nan4, nan5, nan6, neginf, posinf, negzro, 
    newzro;

  ret_val = 1;

  posinf = *one / *zero;
  if (posinf <= *one) {
    ret_val = 0;
    return ret_val;
  }

  neginf = -(*one) / *zero;
  if (neginf >= *zero) {
    ret_val = 0;
    return ret_val;
  }

  negzro = *one / (neginf + *one);
  if (negzro != *zero) {
    ret_val = 0;
    return ret_val;
  }

  neginf = *one / negzro;
  if (neginf >= *zero) {
    ret_val = 0;
    return ret_val;
  }

  newzro = negzro + *zero;
  if (newzro != *zero) {
    ret_val = 0;
    return ret_val;
  }

  posinf = *one / newzro;
  if (posinf <= *one) {
    ret_val = 0;
    return ret_val;
  }

  neginf *= posinf;
  if (neginf >= *zero) {
    ret_val = 0;
    return ret_val;
  }

  posinf *= posinf;
  if (posinf <= *one) {
    ret_val = 0;
    return ret_val;
  }

  /*     Return if we were only asked to check infinity arithmetic */
  if (*ispec == 0) {
    return ret_val;
  }

  nan1 = posinf + neginf;

  nan2 = posinf / neginf;

  nan3 = posinf / posinf;

  nan4 = posinf * *zero;

  nan5 = neginf * negzro;

  nan6 = nan5 * *zero;

  if (nan1 == nan1) {
    ret_val = 0;
    return ret_val;
  }

  if (nan2 == nan2) {
    ret_val = 0;
    return ret_val;
  }

  if (nan3 == nan3) {
    ret_val = 0;
    return ret_val;
  }

  if (nan4 == nan4) {
    ret_val = 0;
    return ret_val;
  }

  if (nan5 == nan5) {
    ret_val = 0;
    return ret_val;
  }

  if (nan6 == nan6) {
    ret_val = 0;
    return ret_val;
  }

  return ret_val;
} /* ieeeck_ */

//-------------------------------------------------------------------------------------------------







// LAPACK auxiliary routine (version 3.8.0) 
integer ilaenv(integer *ispec, char *name__, char *opts, integer *n1, 
  integer *n2, integer *n3, integer *n4, ftnlen name_len, ftnlen opts_len)
{
  /* Table of constant values */

  static integer c__1 = 1;
  static f2c_real c_b173 = 0.f;
  static f2c_real c_b174 = 1.f;
  static integer c__0 = 0;


  /* System generated locals */
  integer ret_val;

  /* Builtin functions */
  /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
  integer i_len(char *, ftnlen), s_cmp(char *, char *, ftnlen, ftnlen);

  /* Local variables */
  static logical twostage;
  static integer i__;
  static char c1[1], c2[2], c3[3], c4[2];
  static integer ic, nb, iz, nx;
  static logical cname;
  static integer nbmin;
  static logical sname;
  extern integer ieeeck(integer *, f2c_real *, f2c_real *);
  static char subnam[16];
  extern integer iparmq(integer *, char *, char *, integer *, integer *, 
    integer *, integer *, ftnlen, ftnlen);


  switch (*ispec) {
  case 1:  goto L10;
  case 2:  goto L10;
  case 3:  goto L10;
  case 4:  goto L80;
  case 5:  goto L90;
  case 6:  goto L100;
  case 7:  goto L110;
  case 8:  goto L120;
  case 9:  goto L130;
  case 10:  goto L140;
  case 11:  goto L150;
  case 12:  goto L160;
  case 13:  goto L160;
  case 14:  goto L160;
  case 15:  goto L160;
  case 16:  goto L160;
  }

  /*     Invalid value for ISPEC */

  ret_val = -1;
  return ret_val;

L10:

  /*     Convert NAME to upper case if the first character is lower case. */

  ret_val = 1;
  s_copy(subnam, name__, (ftnlen)16, name_len);
  ic = *(unsigned char *)subnam;
  iz = 'Z';
  if (iz == 90 || iz == 122) {

    /*        ASCII character set */

    if (ic >= 97 && ic <= 122) {
      *(unsigned char *)subnam = (char) (ic - 32);
      for (i__ = 2; i__ <= 6; ++i__) {
        ic = *(unsigned char *)&subnam[i__ - 1];
        if (ic >= 97 && ic <= 122) {
          *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
        }
        /* L20: */
      }
    }

  } else if (iz == 233 || iz == 169) {

    /*        EBCDIC character set */

    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 && 
      ic <= 169) {
      *(unsigned char *)subnam = (char) (ic + 64);
      for (i__ = 2; i__ <= 6; ++i__) {
        ic = *(unsigned char *)&subnam[i__ - 1];
        if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 
          162 && ic <= 169) {
          *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
        }
        /* L30: */
      }
    }

  } else if (iz == 218 || iz == 250) {

    /*        Prime machines:  ASCII+128 */

    if (ic >= 225 && ic <= 250) {
      *(unsigned char *)subnam = (char) (ic - 32);
      for (i__ = 2; i__ <= 6; ++i__) {
        ic = *(unsigned char *)&subnam[i__ - 1];
        if (ic >= 225 && ic <= 250) {
          *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
        }
        /* L40: */
      }
    }
  }

  *(unsigned char *)c1 = *(unsigned char *)subnam;
  sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
  cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
  if (! (cname || sname)) {
    return ret_val;
  }
  s_copy(c2, subnam + 1, (ftnlen)2, (ftnlen)2);
  s_copy(c3, subnam + 3, (ftnlen)3, (ftnlen)3);
  s_copy(c4, c3 + 1, (ftnlen)2, (ftnlen)2);
  twostage = i_len(subnam, (ftnlen)16) >= 11 && *(unsigned char *)&subnam[
    10] == '2';

  switch (*ispec) {
  case 1:  goto L50;
  case 2:  goto L60;
  case 3:  goto L70;
  }

L50:

  /*     ISPEC = 1:  block size */

  /*     In these examples, separate code is provided for setting NB for */
  /*     real and complex.  We assume that NB will take the same value in */
  /*     single or double precision. */

  nb = 1;

  if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    } else if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, 
      "RQF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)
        3, (ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) 
      == 0) {
      if (sname) {
        nb = 32;
      } else {
        nb = 32;
      }
    } else if (s_cmp(c3, "QR ", (ftnlen)3, (ftnlen)3) == 0) {
      if (*n3 == 1) {
        if (sname) {
          /*     M*N */
          if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
            nb = *n1;
          } else {
            nb = 32768 / *n2;
          }
        } else {
          if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
            nb = *n1;
          } else {
            nb = 32768 / *n2;
          }
        }
      } else {
        if (sname) {
          nb = 1;
        } else {
          nb = 1;
        }
      }
    } else if (s_cmp(c3, "LQ ", (ftnlen)3, (ftnlen)3) == 0) {
      if (*n3 == 2) {
        if (sname) {
          /*     M*N */
          if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
            nb = *n1;
          } else {
            nb = 32768 / *n2;
          }
        } else {
          if (*n1 * *n2 <= 131072 || *n1 <= 8192) {
            nb = *n1;
          } else {
            nb = 32768 / *n2;
          }
        }
      } else {
        if (sname) {
          nb = 1;
        } else {
          nb = 1;
        }
      }
    } else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 32;
      } else {
        nb = 32;
      }
    } else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 32;
      } else {
        nb = 32;
      }
    } else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    }
  } else if (s_cmp(c2, "PO", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    }
  } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        if (twostage) {
          nb = 192;
        } else {
          nb = 64;
        }
      } else {
        if (twostage) {
          nb = 192;
        } else {
          nb = 64;
        }
      }
    } else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nb = 32;
    } else if (sname && s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
      nb = 64;
    }
  } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (twostage) {
        nb = 192;
      } else {
        nb = 64;
      }
    } else if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nb = 32;
    } else if (s_cmp(c3, "GST", (ftnlen)3, (ftnlen)3) == 0) {
      nb = 64;
    }
  } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nb = 32;
      }
    } else if (*(unsigned char *)c3 == 'M') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nb = 32;
      }
    }
  } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nb = 32;
      }
    } else if (*(unsigned char *)c3 == 'M') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nb = 32;
      }
    }
  } else if (s_cmp(c2, "GB", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        if (*n4 <= 64) {
          nb = 1;
        } else {
          nb = 32;
        }
      } else {
        if (*n4 <= 64) {
          nb = 1;
        } else {
          nb = 32;
        }
      }
    }
  } else if (s_cmp(c2, "PB", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        if (*n2 <= 64) {
          nb = 1;
        } else {
          nb = 32;
        }
      } else {
        if (*n2 <= 64) {
          nb = 1;
        } else {
          nb = 32;
        }
      }
    }
  } else if (s_cmp(c2, "TR", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    } else if (s_cmp(c3, "EVC", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    }
  } else if (s_cmp(c2, "LA", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "UUM", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 64;
      } else {
        nb = 64;
      }
    }
  } else if (sname && s_cmp(c2, "ST", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "EBZ", (ftnlen)3, (ftnlen)3) == 0) {
      nb = 1;
    }
  } else if (s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0) {
    nb = 32;
    if (s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nb = 32;
      } else {
        nb = 32;
      }
    }
  }
  ret_val = nb;
  return ret_val;

L60:

  /*     ISPEC = 2:  minimum block size */

  nbmin = 2;
  if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
      ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
        ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
    {
      if (sname) {
        nbmin = 2;
      } else {
        nbmin = 2;
      }
    } else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nbmin = 2;
      } else {
        nbmin = 2;
      }
    } else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nbmin = 2;
      } else {
        nbmin = 2;
      }
    } else if (s_cmp(c3, "TRI", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nbmin = 2;
      } else {
        nbmin = 2;
      }
    }
  } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRF", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nbmin = 8;
      } else {
        nbmin = 8;
      }
    } else if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nbmin = 2;
    }
  } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nbmin = 2;
    }
  } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nbmin = 2;
      }
    } else if (*(unsigned char *)c3 == 'M') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nbmin = 2;
      }
    }
  } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nbmin = 2;
      }
    } else if (*(unsigned char *)c3 == 'M') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nbmin = 2;
      }
    }
  } else if (s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0) {
    nbmin = 2;
    if (s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0) {
      nbmin = 2;
    }
  }
  ret_val = nbmin;
  return ret_val;

L70:

  /*     ISPEC = 3:  crossover point */

  nx = 0;
  if (s_cmp(c2, "GE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "QRF", (ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "RQF", (
      ftnlen)3, (ftnlen)3) == 0 || s_cmp(c3, "LQF", (ftnlen)3, (
        ftnlen)3) == 0 || s_cmp(c3, "QLF", (ftnlen)3, (ftnlen)3) == 0)
    {
      if (sname) {
        nx = 128;
      } else {
        nx = 128;
      }
    } else if (s_cmp(c3, "HRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nx = 128;
      } else {
        nx = 128;
      }
    } else if (s_cmp(c3, "BRD", (ftnlen)3, (ftnlen)3) == 0) {
      if (sname) {
        nx = 128;
      } else {
        nx = 128;
      }
    }
  } else if (s_cmp(c2, "SY", (ftnlen)2, (ftnlen)2) == 0) {
    if (sname && s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nx = 32;
    }
  } else if (cname && s_cmp(c2, "HE", (ftnlen)2, (ftnlen)2) == 0) {
    if (s_cmp(c3, "TRD", (ftnlen)3, (ftnlen)3) == 0) {
      nx = 32;
    }
  } else if (sname && s_cmp(c2, "OR", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nx = 128;
      }
    }
  } else if (cname && s_cmp(c2, "UN", (ftnlen)2, (ftnlen)2) == 0) {
    if (*(unsigned char *)c3 == 'G') {
      if (s_cmp(c4, "QR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "RQ", 
        (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "LQ", (ftnlen)2, (
          ftnlen)2) == 0 || s_cmp(c4, "QL", (ftnlen)2, (ftnlen)2) ==
        0 || s_cmp(c4, "HR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(
          c4, "TR", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c4, "BR", (
            ftnlen)2, (ftnlen)2) == 0) {
        nx = 128;
      }
    }
  } else if (s_cmp(c2, "GG", (ftnlen)2, (ftnlen)2) == 0) {
    nx = 128;
    if (s_cmp(c3, "HD3", (ftnlen)3, (ftnlen)3) == 0) {
      nx = 128;
    }
  }
  ret_val = nx;
  return ret_val;

L80:

  /*     ISPEC = 4:  number of shifts (used by xHSEQR) */

  ret_val = 6;
  return ret_val;

L90:

  /*     ISPEC = 5:  minimum column dimension (not used) */

  ret_val = 2;
  return ret_val;

L100:

  /*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */

  ret_val = (integer) ((f2c_real) min(*n1,*n2) * 1.6f);
  return ret_val;

L110:

  /*     ISPEC = 7:  number of processors (not used) */

  ret_val = 1;
  return ret_val;

L120:

  /*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

  ret_val = 50;
  return ret_val;

L130:

  /*     ISPEC = 9:  maximum size of the subproblems at the bottom of the */
  /*                 computation tree in the divide-and-conquer algorithm */
  /*                 (used by xGELSD and xGESDD) */

  ret_val = 25;
  return ret_val;

L140:

  /*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap */

  /*     ILAENV = 0 */
  ret_val = 1;
  if (ret_val == 1) {
    ret_val = ieeeck(&c__1, &c_b173, &c_b174);
  }
  return ret_val;

L150:

  /*     ISPEC = 11: infinity arithmetic can be trusted not to trap */

  /*     ILAENV = 0 */
  ret_val = 1;
  if (ret_val == 1) {
    ret_val = ieeeck(&c__0, &c_b173, &c_b174);
  }
  return ret_val;

L160:

  /*     12 <= ISPEC <= 16: xHSEQR or related subroutines. */

  ret_val = iparmq(ispec, name__, opts, n1, n2, n3, n4, name_len, opts_len)
    ;
  return ret_val;

  /*     End of ILAENV */

} /* ilaenv_ */

//-------------------------------------------------------------------------------------------------



// LAPACK auxiliary routine (version 3.7.1)
integer iparmq(integer *ispec, char *name__, char *opts, integer *n, integer 
  *ilo, integer *ihi, integer *lwork, ftnlen name_len, ftnlen opts_len)
{
  /* System generated locals */
  integer ret_val, i__1, i__2;
  f2c_real r__1;

  /* Builtin functions */
  double log(doublereal);
  integer i_nint(f2c_real *);
  /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
  integer s_cmp(char *, char *, ftnlen, ftnlen);

  /* Local variables */
  static integer i__, ic, nh, ns, iz;
  static char subnam[6];

  if (*ispec == 15 || *ispec == 13 || *ispec == 16) {

    /*        ==== Set the number simultaneous shifts ==== */

    nh = *ihi - *ilo + 1;
    ns = 2;
    if (nh >= 30) {
      ns = 4;
    }
    if (nh >= 60) {
      ns = 10;
    }
    if (nh >= 150) {
      /* Computing MAX */
      r__1 = (f2c_real) (log((f2c_real) nh) / log(2.f)); // outer cast to f2c_real added by Robin Schmidt
      i__1 = 10, i__2 = nh / i_nint(&r__1);
      ns = max(i__1,i__2);
    }
    if (nh >= 590) {
      ns = 64;
    }
    if (nh >= 3000) {
      ns = 128;
    }
    if (nh >= 6000) {
      ns = 256;
    }
    /* Computing MAX */
    i__1 = 2, i__2 = ns - ns % 2;
    ns = max(i__1,i__2);
  }

  if (*ispec == 12) {


    /*        ===== Matrices of order smaller than NMIN get sent */
    /*        .     to xLAHQR, the classic double shift algorithm. */
    /*        .     This must be at least 11. ==== */

    ret_val = 75;

  } else if (*ispec == 14) {

    /*        ==== INIBL: skip a multi-shift qr iteration and */
    /*        .    whenever aggressive early deflation finds */
    /*        .    at least (NIBBLE*(window size)/100) deflations. ==== */

    ret_val = 14;

  } else if (*ispec == 15) {

    /*        ==== NSHFTS: The number of simultaneous shifts ===== */

    ret_val = ns;

  } else if (*ispec == 13) {

    /*        ==== NW: deflation window size.  ==== */

    if (nh <= 500) {
      ret_val = ns;
    } else {
      ret_val = ns * 3 / 2;
    }

  } else if (*ispec == 16) {

    /*        ==== IACC22: Whether to accumulate reflections */
    /*        .     before updating the far-from-diagonal elements */
    /*        .     and whether to use 2-by-2 block structure while */
    /*        .     doing it.  A small amount of work could be saved */
    /*        .     by making this choice dependent also upon the */
    /*        .     NH=IHI-ILO+1. */


    /*        Convert NAME to upper case if the first character is lower case. */

    ret_val = 0;
    s_copy(subnam, name__, (ftnlen)6, name_len);
    ic = *(unsigned char *)subnam;
    iz = 'Z';
    if (iz == 90 || iz == 122) {

      /*           ASCII character set */

      if (ic >= 97 && ic <= 122) {
        *(unsigned char *)subnam = (char) (ic - 32);
        for (i__ = 2; i__ <= 6; ++i__) {
          ic = *(unsigned char *)&subnam[i__ - 1];
          if (ic >= 97 && ic <= 122) {
            *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
          }
        }
      }

    } else if (iz == 233 || iz == 169) {

      /*           EBCDIC character set */

      if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 
        && ic <= 169) {
        *(unsigned char *)subnam = (char) (ic + 64);
        for (i__ = 2; i__ <= 6; ++i__) {
          ic = *(unsigned char *)&subnam[i__ - 1];
          if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || 
            ic >= 162 && ic <= 169) {
            *(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
          }
        }
      }

    } else if (iz == 218 || iz == 250) {

      /*           Prime machines:  ASCII+128 */

      if (ic >= 225 && ic <= 250) {
        *(unsigned char *)subnam = (char) (ic - 32);
        for (i__ = 2; i__ <= 6; ++i__) {
          ic = *(unsigned char *)&subnam[i__ - 1];
          if (ic >= 225 && ic <= 250) {
            *(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
          }
        }
      }
    }

    if (s_cmp(subnam + 1, "GGHRD", (ftnlen)5, (ftnlen)5) == 0 || s_cmp(
      subnam + 1, "GGHD3", (ftnlen)5, (ftnlen)5) == 0) {
      ret_val = 1;
      if (nh >= 14) {
        ret_val = 2;
      }
    } else if (s_cmp(subnam + 3, "EXC", (ftnlen)3, (ftnlen)3) == 0) {
      if (nh >= 14) {
        ret_val = 1;
      }
      if (nh >= 14) {
        ret_val = 2;
      }
    } else if (s_cmp(subnam + 1, "HSEQR", (ftnlen)5, (ftnlen)5) == 0 || 
      s_cmp(subnam + 1, "LAQR", (ftnlen)4, (ftnlen)4) == 0) {
      if (ns >= 14) {
        ret_val = 1;
      }
      if (ns >= 14) {
        ret_val = 2;
      }
    }

  } else {
    /*        ===== invalid value of ispec ===== */
    ret_val = -1;

  }

  /*     ==== End of IPARMQ ==== */

  return ret_val;
} /* iparmq_ */

//-------------------------------------------------------------------------------------------------

// translated from dlaswp, LAPACK auxiliary routine (version 3.7.1)
template<class T>
int laswp(integer *n, T *a, integer *lda, integer *k1, integer *k2, integer *ipiv, integer *incx)
{
  /* System generated locals */
  integer a_dim1, a_offset, i__1, i__2, i__3, i__4;

  /* Local variables */
  static integer i__, j, k, i1, i2, n32, ip, ix, ix0, inc;
  static T temp;

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  --ipiv;

  /* Function Body */
  if (*incx > 0) {
    ix0 = *k1;
    i1 = *k1;
    i2 = *k2;
    inc = 1;
  } else if (*incx < 0) {
    ix0 = *k1 + (*k1 - *k2) * *incx;
    i1 = *k2;
    i2 = *k1;
    inc = -1;
  } else {
    return 0;
  }

  n32 = *n / 32 << 5;
  if (n32 != 0) {
    i__1 = n32;
    for (j = 1; j <= i__1; j += 32) {
      ix = ix0;
      i__2 = i2;
      i__3 = inc;
      for (i__ = i1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) 
      {
        ip = ipiv[ix];
        if (ip != i__) {
          i__4 = j + 31;
          for (k = j; k <= i__4; ++k) {
            temp = a[i__ + k * a_dim1];
            a[i__ + k * a_dim1] = a[ip + k * a_dim1];
            a[ip + k * a_dim1] = temp;
            /* L10: */
          }
        }
        ix += *incx;
        /* L20: */
      }
      /* L30: */
    }
  }
  if (n32 != *n) {
    ++n32;
    ix = ix0;
    i__1 = i2;
    i__3 = inc;
    for (i__ = i1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
      ip = ipiv[ix];
      if (ip != i__) {
        i__2 = *n;
        for (k = n32; k <= i__2; ++k) {
          temp = a[i__ + k * a_dim1];
          a[i__ + k * a_dim1] = a[ip + k * a_dim1];
          a[ip + k * a_dim1] = temp;
          /* L40: */
        }
      }
      ix += *incx;
      /* L50: */
    }
  }

  return 0;

  /*     End of DLASWP */

} /* dlaswp_ */




// commented because of linker errors - now all these subroutines:
// dgemm*, dcopy*, dswap, dtrsm, idamax, ilaenv, dgemv, dtbsv, dgbtf2, dlaswp,  dscal, 
// have to translated to fix the linker errors

template int gbtrf(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab, 
  integer *ipiv, integer *info);

template int gbtrs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, double *ab,
    integer *ldab, integer *ipiv, double *b, integer *ldb, integer *info, ftnlen trans_len);

/*
template int gbsv(long int *n, long int *kl, long int *ku, long int *nrhs, double *ab, 
  long int *ldab, long int *ipiv, double *b, long int *ldb, long int *info);
*/








}