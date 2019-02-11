#include "LaPack.hpp"
namespace LaPackCPP {

//  -- LAPACK driver routine (version 3.7.0) -- 
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

//  -- LAPACK computational routine (version 3.7.0) -- 
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
  extern int dger_(integer *, integer *, T *,
    T *, integer *, T *, integer *, T *,
    integer *);
  static T temp;
  extern int dscal_(integer *, T *, T *,
    integer *), dgemm_(char *, char *, integer *, integer *, integer *
      , T *, T *, integer *, T *, integer *,
      T *, T *, integer *, ftnlen, ftnlen), dcopy_(
        integer *, T *, integer *, T *, integer *),
    dswap_(integer *, T *, integer *, T *, integer *
      );
  static T work13[4160], work31[4160];
  extern int dtrsm_(char *, char *, char *, char *,
    integer *, integer *, T *, T *, integer *,
    T *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dgbtf2_(
      integer *, integer *, integer *, integer *, T *, integer
      *, integer *, integer *);
  extern integer idamax_(integer *, T *, integer *);
  extern int xerbla_(char *, integer *, ftnlen);
  extern integer ilaenv_(integer *, char *, char *, integer *, integer *,
    integer *, integer *, ftnlen, ftnlen);
  extern int dlaswp_(integer *, T *, integer *,
    integer *, integer *, integer *, integer *);


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
    xerbla_("DGBTRF", &i__1, (ftnlen)6);
    return 0;
  }

  // Quick return if possible
  if(*m == 0 || *n == 0) {
    return 0;
  }

  // Determine the block size for this environment
  nb = ilaenv_(&c__1, "DGBTRF", " ", m, n, kl, ku, (ftnlen)6, (ftnlen)1);
  // The block size must not exceed the limit set by the size of the
  // local arrays WORK13 and WORK31.

  nb = min(nb, 64);

  if(nb <= 1 || nb > *kl) {
    // Use unblocked code
    dgbtf2_(m, n, kl, ku, &ab[ab_offset], ldab, &ipiv[1], info);
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
        jp = idamax_(&i__4, &ab[kv + 1 + jj * ab_dim1], &c__1);
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
              dswap_(&jb, &ab[kv + 1 + jj - j + j * ab_dim1], &
                i__4, &ab[kv + jp + jj - j + j * ab_dim1],
                &i__5);
            }
            else {

            // The interchange affects columns J to JJ-1 of A31 
            // which are stored in the work array WORK31 
              i__4 = jj - j;
              i__5 = *ldab - 1;
              dswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1],
                &i__5, &work31[jp + jj - j - *kl - 1], &
                c__65);
              i__4 = j + jb - jj;
              i__5 = *ldab - 1;
              i__6 = *ldab - 1;
              dswap_(&i__4, &ab[kv + 1 + jj * ab_dim1], &i__5, &
                ab[kv + jp + jj * ab_dim1], &i__6);
            }
          }

          // Compute multipliers 

          d__1 = 1. / ab[kv + 1 + jj * ab_dim1];
          dscal_(&km, &d__1, &ab[kv + 2 + jj * ab_dim1], &c__1);

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
            dger_(&km, &i__4, &c_b18, &ab[kv + 2 + jj * ab_dim1],
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
          dcopy_(&nw, &ab[kv + *kl + 1 - jj + j + jj * ab_dim1], &
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
        dlaswp_(&j2, &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__3, &
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
          dtrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j2,
            &c_b31, &ab[kv + 1 + j * ab_dim1], &i__3, &ab[kv
            + 1 - jb + (j + jb) * ab_dim1], &i__4, (ftnlen)4,
            (ftnlen)5, (ftnlen)12, (ftnlen)4);

          if(i2 > 0) {

            // Update A22
            i__3 = *ldab - 1;
            i__4 = *ldab - 1;
            i__5 = *ldab - 1;
            dgemm_("No transpose", "No transpose", &i2, &j2, &jb,
              &c_b18, &ab[kv + 1 + jb + j * ab_dim1], &i__3,
              &ab[kv + 1 - jb + (j + jb) * ab_dim1], &i__4,
              &c_b31, &ab[kv + 1 + (j + jb) * ab_dim1], &
              i__5, (ftnlen)12, (ftnlen)12);
          }

          if(i3 > 0) {

            // Update A32
            i__3 = *ldab - 1;
            i__4 = *ldab - 1;
            dgemm_("No transpose", "No transpose", &i3, &j2, &jb,
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
          dtrsm_("Left", "Lower", "No transpose", "Unit", &jb, &j3,
            &c_b31, &ab[kv + 1 + j * ab_dim1], &i__3, work13,
            &c__65, (ftnlen)4, (ftnlen)5, (ftnlen)12, (ftnlen)
            4);

          if(i2 > 0) {
            // Update A23 
            i__3 = *ldab - 1;
            i__4 = *ldab - 1;
            dgemm_("No transpose", "No transpose", &i2, &j3, &jb,
              &c_b18, &ab[kv + 1 + jb + j * ab_dim1], &i__3,
              work13, &c__65, &c_b31, &ab[jb + 1 + (j + kv)
              * ab_dim1], &i__4, (ftnlen)12, (ftnlen)12);
          }

          if(i3 > 0) {
            // Update A33 
            i__3 = *ldab - 1;
            dgemm_("No transpose", "No transpose", &i3, &j3, &jb,
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
            dswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
              i__5, &ab[kv + jp + jj - j + j * ab_dim1], &
              i__6);
          }
          else {
            // The interchange does affect A31
            i__4 = jj - j;
            i__5 = *ldab - 1;
            dswap_(&i__4, &ab[kv + 1 + jj - j + j * ab_dim1], &
              i__5, &work31[jp + jj - j - *kl - 1], &c__65);
          }
        }

        // Copy the current column of A31 back into place
        // Computing MIN 
        i__4 = i3, i__5 = jj - j + 1;
        nw = min(i__4, i__5);
        if(nw > 0) {
          dcopy_(&nw, &work31[(jj - j + 1) * 65 - 65], &c__1, &ab[
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

// -- LAPACK computational routine (version 3.7.0) --
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
  extern int dger_(integer *, integer *, doublereal *,
    doublereal *, integer *, doublereal *, integer *, doublereal *,
    integer *);
  extern logical lsame_(char *, char *, ftnlen, ftnlen);
  extern int dgemv_(char *, integer *, integer *,
    doublereal *, doublereal *, integer *, doublereal *, integer *,
    doublereal *, doublereal *, integer *, ftnlen), dswap_(integer *,
      doublereal *, integer *, doublereal *, integer *), dtbsv_(char *,
        char *, char *, integer *, integer *, doublereal *, integer *,
        doublereal *, integer *, ftnlen, ftnlen, ftnlen);
  static logical lnoti;
  extern int xerbla(char *, integer *, ftnlen);
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
  notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
  if(!notran && !lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && !lsame_(
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
          dswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
        }
        dger_(&lm, nrhs, &c_b7, &ab[kd + 1 + j * ab_dim1], &c__1, &b[
          j + b_dim1], ldb, &b[j + 1 + b_dim1], ldb);
        // L10:
      }
    }

    i__1 = *nrhs;
    for(i__ = 1; i__ <= i__1; ++i__) {

      // Solve U*X = B, overwriting B with X.
      i__2 = *kl + *ku;
      dtbsv_("Upper", "No transpose", "Non-unit", n, &i__2, &ab[
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
      dtbsv_("Upper", "Transpose", "Non-unit", n, &i__2, &ab[ab_offset],
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
        dgemv_("Transpose", &lm, nrhs, &c_b7, &b[j + 1 + b_dim1], ldb,
          &ab[kd + 1 + j * ab_dim1], &c__1, &c_b23, &b[j +
          b_dim1], ldb, (ftnlen)9);
        l = ipiv[j];
        if(l != j) {
          dswap_(nrhs, &b[l + b_dim1], ldb, &b[j + b_dim1], ldb);
        }
        // L40:
      }
    }
  }
  return 0;

  // End of DGBTRS 

} // gbtrs




// commented because of linker errors - now all these subroutines:
// dger, dgemm, dcopy, dswap, dtrsm, idamax, ilaenv, dgemv, dtbsv, dgbtf2, dlaswp,  dscal, 
// have to translated to fix the linker errors

/*
template int gbtrf(integer *m, integer *n, integer *kl, integer *ku, double *ab, integer *ldab, 
  integer *ipiv, integer *info);

template int gbtrs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, double *ab,
    integer *ldab, integer *ipiv, double *b, integer *ldb, integer *info, ftnlen trans_len);

template int gbsv(long int *n, long int *kl, long int *ku, long int *nrhs, double *ab, 
  long int *ldab, long int *ipiv, double *b, long int *ldb, long int *info);
*/








}