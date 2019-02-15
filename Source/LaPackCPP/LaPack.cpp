
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
// DRIVER routines

// translated from dgbsv, LAPACK driver routine (version 3.7.0) -- 
template<class T>
int gbsv(long int *n, long int *kl, long int *ku, long int *nrhs, T *ab, long int *ldab,
  long int *ipiv, T *b, long int *ldb, long int *info)
{
  // System generated locals
  long int ab_dim1, ab_offset, b_dim1, b_offset, i__1;

  // Local variables
  // Subroutine
  //extern  int
  //  gbtrf(long int *, long int *, long int *,
  //    long int *, double *, long int *, long int *, long int *),
  //  xerbla_(char *, long int *, ftnlen),
  //  gbtrs(char *, long int *,
  //    long int *, long int *, long int *, double *, long int *, long int *, double *,
  //    long int *, long int *, ftnlen);

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
    xerbla("DGBSV ", &i__1, (ftnlen)6);
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

// translated from dgbsvx - LAPACK driver routine (version 3.7.0)
template<class T>
int gbsvx(char *fact, char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, 
  T *ab, integer *ldab, T *afb, integer *ldafb, integer *ipiv, char *equed, T *r__, T *c__, T *b, 
  integer *ldb, T *x, integer *ldx, T *rcond, T *ferr, T *berr, T *work, integer *iwork, 
  integer *info, ftnlen fact_len, ftnlen trans_len, ftnlen equed_len)
{
  /* System generated locals */
  integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5;
  T d__1, d__2, d__3;

  /* Local variables */
  static integer i__, j, j1, j2;
  static T amax;
  static char norm[1];
  //extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static T rcmin, rcmax, anorm;
  //extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
  //  doublereal *, integer *);
  static logical equil;
  extern doublereal dlangb_(char *, integer *, integer *, integer *, 
    doublereal *, integer *, doublereal *, ftnlen), dlamch_(char *, 
      ftnlen);
  extern /* Subroutine */ int dlaqgb_(integer *, integer *, integer *, 
    integer *, doublereal *, integer *, doublereal *, doublereal *, 
    doublereal *, doublereal *, doublereal *, char *, ftnlen), 
    dgbcon_(char *, integer *, integer *, integer *, doublereal *, 
      integer *, integer *, doublereal *, doublereal *, doublereal *, 
      integer *, integer *, ftnlen);
  static T colcnd;
  extern doublereal dlantb_(char *, char *, char *, integer *, integer *, 
    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
  extern /* Subroutine */ int dgbequ_(integer *, integer *, integer *, 
    integer *, doublereal *, integer *, doublereal *, doublereal *, 
    doublereal *, doublereal *, doublereal *, integer *), dgbrfs_(
      char *, integer *, integer *, integer *, integer *, doublereal *, 
      integer *, doublereal *, integer *, integer *, doublereal *, 
      integer *, doublereal *, integer *, doublereal *, doublereal *, 
      doublereal *, integer *, integer *, ftnlen), dgbtrf_(integer *, 
        integer *, integer *, integer *, doublereal *, integer *, integer 
        *, integer *);
  static logical nofact;
  extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
    doublereal *, integer *, doublereal *, integer *, ftnlen), 
    xerbla_(char *, integer *, ftnlen);
  static T bignum;
  extern /* Subroutine */ int dgbtrs_(char *, integer *, integer *, integer 
    *, integer *, doublereal *, integer *, integer *, doublereal *, 
    integer *, integer *, ftnlen);
  static integer infequ;
  static logical colequ;
  static T rowcnd;
  static logical notran;
  static T smlnum;
  static logical rowequ;
  static T rpvgrw;/

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  afb_dim1 = *ldafb;
  afb_offset = 1 + afb_dim1;
  afb -= afb_offset;
  --ipiv;
  --r__;
  --c__;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;
  x_dim1 = *ldx;
  x_offset = 1 + x_dim1;
  x -= x_offset;
  --ferr;
  --berr;
  --work;
  --iwork;

  /* Function Body */
  *info = 0;
  nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
  equil = lsame_(fact, "E", (ftnlen)1, (ftnlen)1);
  notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
  if (nofact || equil) {
    *(unsigned char *)equed = 'N';
    rowequ = FALSE_;
    colequ = FALSE_;
  } else {
    rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
      "B", (ftnlen)1, (ftnlen)1);
    colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed, 
      "B", (ftnlen)1, (ftnlen)1);
    smlnum = dlamch_("Safe minimum", (ftnlen)12);
    bignum = 1. / smlnum;
  }

  /*     Test the input parameters. */

  if (! nofact && ! equil && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
    *info = -1;
  } else if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! 
    lsame_(trans, "C", (ftnlen)1, (ftnlen)1)) {
    *info = -2;
  } else if (*n < 0) {
    *info = -3;
  } else if (*kl < 0) {
    *info = -4;
  } else if (*ku < 0) {
    *info = -5;
  } else if (*nrhs < 0) {
    *info = -6;
  } else if (*ldab < *kl + *ku + 1) {
    *info = -8;
  } else if (*ldafb < (*kl << 1) + *ku + 1) {
    *info = -10;
  } else if (lsame_(fact, "F", (ftnlen)1, (ftnlen)1) && ! (rowequ || colequ 
    || lsame_(equed, "N", (ftnlen)1, (ftnlen)1))) {
    *info = -12;
  } else {
    if (rowequ) {
      rcmin = bignum;
      rcmax = 0.;
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        /* Computing MIN */
        d__1 = rcmin, d__2 = r__[j];
        rcmin = min(d__1,d__2);
        /* Computing MAX */
        d__1 = rcmax, d__2 = r__[j];
        rcmax = max(d__1,d__2);
        /* L10: */
      }
      if (rcmin <= 0.) {
        *info = -13;
      } else if (*n > 0) {
        rowcnd = max(rcmin,smlnum) / min(rcmax,bignum);
      } else {
        rowcnd = 1.;
      }
    }
    if (colequ && *info == 0) {
      rcmin = bignum;
      rcmax = 0.;
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        /* Computing MIN */
        d__1 = rcmin, d__2 = c__[j];
        rcmin = min(d__1,d__2);
        /* Computing MAX */
        d__1 = rcmax, d__2 = c__[j];
        rcmax = max(d__1,d__2);
        /* L20: */
      }
      if (rcmin <= 0.) {
        *info = -14;
      } else if (*n > 0) {
        colcnd = max(rcmin,smlnum) / min(rcmax,bignum);
      } else {
        colcnd = 1.;
      }
    }
    if (*info == 0) {
      if (*ldb < max(1,*n)) {
        *info = -16;
      } else if (*ldx < max(1,*n)) {
        *info = -18;
      }
    }
  }

  if (*info != 0) {
    i__1 = -(*info);
    xerbla_("DGBSVX", &i__1, (ftnlen)6);
    return 0;
  }

  if (equil) {

    /*        Compute row and column scalings to equilibrate the matrix A. */

    dgbequ_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &rowcnd,
      &colcnd, &amax, &infequ);
    if (infequ == 0) {

      /*           Equilibrate the matrix. */

      dlaqgb_(n, n, kl, ku, &ab[ab_offset], ldab, &r__[1], &c__[1], &
        rowcnd, &colcnd, &amax, equed, (ftnlen)1);
      rowequ = lsame_(equed, "R", (ftnlen)1, (ftnlen)1) || lsame_(equed,
        "B", (ftnlen)1, (ftnlen)1);
      colequ = lsame_(equed, "C", (ftnlen)1, (ftnlen)1) || lsame_(equed,
        "B", (ftnlen)1, (ftnlen)1);
    }
  }

  /*     Scale the right hand side. */

  if (notran) {
    if (rowequ) {
      i__1 = *nrhs;
      for (j = 1; j <= i__1; ++j) {
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          b[i__ + j * b_dim1] = r__[i__] * b[i__ + j * b_dim1];
          /* L30: */
        }
        /* L40: */
      }
    }
  } else if (colequ) {
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__) {
        b[i__ + j * b_dim1] = c__[i__] * b[i__ + j * b_dim1];
        /* L50: */
      }
      /* L60: */
    }
  }

  if (nofact || equil) {

    /*        Compute the LU factorization of the band matrix A. */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      /* Computing MAX */
      i__2 = j - *ku;
      j1 = max(i__2,1);
      /* Computing MIN */
      i__2 = j + *kl;
      j2 = min(i__2,*n);
      i__2 = j2 - j1 + 1;
      dcopy_(&i__2, &ab[*ku + 1 - j + j1 + j * ab_dim1], &c__1, &afb[*
        kl + *ku + 1 - j + j1 + j * afb_dim1], &c__1);
      /* L70: */
    }

    dgbtrf_(n, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], info);

    /*        Return if INFO is non-zero. */

    if (*info > 0) {

      /*           Compute the reciprocal pivot growth factor of the */
      /*           leading rank-deficient INFO columns of A. */

      anorm = 0.;
      i__1 = *info;
      for (j = 1; j <= i__1; ++j) {
        /* Computing MAX */
        i__2 = *ku + 2 - j;
        /* Computing MIN */
        i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
        i__3 = min(i__4,i__5);
        for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
          /* Computing MAX */
          d__2 = anorm, d__3 = (d__1 = ab[i__ + j * ab_dim1], abs(
            d__1));
          anorm = max(d__2,d__3);
          /* L80: */
        }
        /* L90: */
      }
      /* Computing MIN */
      i__3 = *info - 1, i__2 = *kl + *ku;
      i__1 = min(i__3,i__2);
      /* Computing MAX */
      i__4 = 1, i__5 = *kl + *ku + 2 - *info;
      rpvgrw = dlantb_("M", "U", "N", info, &i__1, &afb[max(i__4,i__5) 
        + afb_dim1], ldafb, &work[1], (ftnlen)1, (ftnlen)1, (
          ftnlen)1);
      if (rpvgrw == 0.) {
        rpvgrw = 1.;
      } else {
        rpvgrw = anorm / rpvgrw;
      }
      work[1] = rpvgrw;
      *rcond = 0.;
      return 0;
    }
  }

  /*     Compute the norm of the matrix A and the */
  /*     reciprocal pivot growth factor RPVGRW. */

  if (notran) {
    *(unsigned char *)norm = '1';
  } else {
    *(unsigned char *)norm = 'I';
  }
  anorm = dlangb_(norm, n, kl, ku, &ab[ab_offset], ldab, &work[1], (ftnlen)
    1);
  i__1 = *kl + *ku;
  rpvgrw = dlantb_("M", "U", "N", n, &i__1, &afb[afb_offset], ldafb, &work[
    1], (ftnlen)1, (ftnlen)1, (ftnlen)1);
  if (rpvgrw == 0.) {
    rpvgrw = 1.;
  } else {
    rpvgrw = dlangb_("M", n, kl, ku, &ab[ab_offset], ldab, &work[1], (
      ftnlen)1) / rpvgrw;
  }

  /*     Compute the reciprocal of the condition number of A. */

  dgbcon_(norm, n, kl, ku, &afb[afb_offset], ldafb, &ipiv[1], &anorm, rcond,
    &work[1], &iwork[1], info, (ftnlen)1);

  /*     Compute the solution matrix X. */

  dlacpy_("Full", n, nrhs, &b[b_offset], ldb, &x[x_offset], ldx, (ftnlen)4);
  dgbtrs_(trans, n, kl, ku, nrhs, &afb[afb_offset], ldafb, &ipiv[1], &x[
    x_offset], ldx, info, (ftnlen)1);

  /*     Use iterative refinement to improve the computed solution and */
  /*     compute error bounds and backward error estimates for it. */

  dgbrfs_(trans, n, kl, ku, nrhs, &ab[ab_offset], ldab, &afb[afb_offset], 
    ldafb, &ipiv[1], &b[b_offset], ldb, &x[x_offset], ldx, &ferr[1], &
    berr[1], &work[1], &iwork[1], info, (ftnlen)1);

  /*     Transform the solution matrix X to a solution of the original */
  /*     system. */

  if (notran) {
    if (colequ) {
      i__1 = *nrhs;
      for (j = 1; j <= i__1; ++j) {
        i__3 = *n;
        for (i__ = 1; i__ <= i__3; ++i__) {
          x[i__ + j * x_dim1] = c__[i__] * x[i__ + j * x_dim1];
          /* L100: */
        }
        /* L110: */
      }
      i__1 = *nrhs;
      for (j = 1; j <= i__1; ++j) {
        ferr[j] /= colcnd;
        /* L120: */
      }
    }
  } else if (rowequ) {
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
      i__3 = *n;
      for (i__ = 1; i__ <= i__3; ++i__) {
        x[i__ + j * x_dim1] = r__[i__] * x[i__ + j * x_dim1];
        /* L130: */
      }
      /* L140: */
    }
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
      ferr[j] /= rowcnd;
      /* L150: */
    }
  }

  /*     Set INFO = N+1 if the matrix is singular to working precision. */

  if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
    *info = *n + 1;
  }

  work[1] = rpvgrw;
  return 0;

} /* gbsvx */


//=================================================================================================
// COMPUTATIONAL routines

// from dgbcon - LAPACK computational routine (version 3.7.0) 
template<class T>
int gbcon(char *norm, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, integer *ipiv, 
  T *anorm, T *rcond, T *work, integer *iwork, integer *info, ftnlen norm_len)
{
  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3;
  T d__1;

  /* Local variables */
  static integer j;
  static T t;
  static integer kd, lm, jp, ix, kase;
  extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
    integer *);
  static integer kase1;
  static T scale;
  extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static integer isave[3];
  extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
    integer *);
  static logical lnoti;
  extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
    integer *, doublereal *, integer *), dlacn2_(integer *, 
      doublereal *, doublereal *, integer *, doublereal *, integer *, 
      integer *);
  extern doublereal dlamch_(char *, ftnlen);
  extern integer idamax_(integer *, doublereal *, integer *);
  extern /* Subroutine */ int dlatbs_(char *, char *, char *, char *, 
    integer *, integer *, doublereal *, integer *, doublereal *, 
    doublereal *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, 
    ftnlen), xerbla_(char *, integer *, ftnlen);
  static T ainvnm;
  static logical onenrm;
  static char normin[1];
  static T smlnum;

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --ipiv;
  --work;
  --iwork;

  /* Function Body */
  *info = 0;
  onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
    ftnlen)1);
  if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
    *info = -1;
  } else if (*n < 0) {
    *info = -2;
  } else if (*kl < 0) {
    *info = -3;
  } else if (*ku < 0) {
    *info = -4;
  } else if (*ldab < (*kl << 1) + *ku + 1) {
    *info = -6;
  } else if (*anorm < 0.) {
    *info = -8;
  }
  if (*info != 0) {
    i__1 = -(*info);
    xerbla_("DGBCON", &i__1, (ftnlen)6);
    return 0;
  }

  /*     Quick return if possible */

  *rcond = 0.;
  if (*n == 0) {
    *rcond = 1.;
    return 0;
  } else if (*anorm == 0.) {
    return 0;
  }

  smlnum = dlamch_("Safe minimum", (ftnlen)12);

  /*     Estimate the norm of inv(A). */

  ainvnm = 0.;
  *(unsigned char *)normin = 'N';
  if (onenrm) {
    kase1 = 1;
  } else {
    kase1 = 2;
  }
  kd = *kl + *ku + 1;
  lnoti = *kl > 0;
  kase = 0;
L10:
  dlacn2_(n, &work[*n + 1], &work[1], &iwork[1], &ainvnm, &kase, isave);
  if (kase != 0) {
    if (kase == kase1) {

      /*           Multiply by inv(L). */

      if (lnoti) {
        i__1 = *n - 1;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MIN */
          i__2 = *kl, i__3 = *n - j;
          lm = min(i__2,i__3);
          jp = ipiv[j];
          t = work[jp];
          if (jp != j) {
            work[jp] = work[j];
            work[j] = t;
          }
          d__1 = -t;
          daxpy_(&lm, &d__1, &ab[kd + 1 + j * ab_dim1], &c__1, &
            work[j + 1], &c__1);
          /* L20: */
        }
      }

      /*           Multiply by inv(U). */

      i__1 = *kl + *ku;
      dlatbs_("Upper", "No transpose", "Non-unit", normin, n, &i__1, &
        ab[ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 
        1], info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
    } else {

      /*           Multiply by inv(U**T). */

      i__1 = *kl + *ku;
      dlatbs_("Upper", "Transpose", "Non-unit", normin, n, &i__1, &ab[
        ab_offset], ldab, &work[1], &scale, &work[(*n << 1) + 1], 
          info, (ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);

      /*           Multiply by inv(L**T). */

      if (lnoti) {
        for (j = *n - 1; j >= 1; --j) {
          /* Computing MIN */
          i__1 = *kl, i__2 = *n - j;
          lm = min(i__1,i__2);
          work[j] -= ddot_(&lm, &ab[kd + 1 + j * ab_dim1], &c__1, &
            work[j + 1], &c__1);
          jp = ipiv[j];
          if (jp != j) {
            t = work[jp];
            work[jp] = work[j];
            work[j] = t;
          }
          /* L30: */
        }
      }
    }

    /*        Divide X by 1/SCALE if doing so will not cause overflow. */

    *(unsigned char *)normin = 'Y';
    if (scale != 1.) {
      ix = idamax_(n, &work[1], &c__1);
      if (scale < (d__1 = work[ix], abs(d__1)) * smlnum || scale == 0.) 
      {
        goto L40;
      }
      drscl_(n, &scale, &work[1], &c__1);
    }
    goto L10;
  }

  /*     Compute the estimate of the reciprocal condition number. */

  if (ainvnm != 0.) {
    *rcond = 1. / ainvnm / *anorm;
  }

L40:
  return 0;

  /*     End of DGBCON */

} /* dgbcon_ */

//-------------------------------------------------------------------------------------------------

// from dgbrfs - LAPACK computational routine (version 3.7.0)
template<class T>
int gbrfs(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, T *ab, integer *ldab,
  T *afb, integer *ldafb, integer *ipiv, T *b, integer *ldb, T *x, integer *ldx, T *ferr, T *berr, 
  T *work, integer *iwork, integer *info, ftnlen trans_len)
{
  /* System generated locals */
  integer ab_dim1, ab_offset, afb_dim1, afb_offset, b_dim1, b_offset, 
    x_dim1, x_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;
  T d__1, d__2, d__3;

  /* Local variables */
  static integer i__, j, k;
  static T s;
  static integer kk;
  static T xk;
  static integer nz;
  static T eps;
  static integer kase;
  static T safe1, safe2;
  extern /* Subroutine */ int dgbmv_(char *, integer *, integer *, integer *
    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
    integer *, doublereal *, doublereal *, integer *, ftnlen);
  extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static integer isave[3];
  extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
    doublereal *, integer *), daxpy_(integer *, doublereal *, 
      doublereal *, integer *, doublereal *, integer *);
  static integer count;
  extern /* Subroutine */ int dlacn2_(integer *, doublereal *, doublereal *,
    integer *, doublereal *, integer *, integer *);
  extern doublereal dlamch_(char *, ftnlen);
  static T safmin;
  extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dgbtrs_(
    char *, integer *, integer *, integer *, integer *, doublereal *, 
    integer *, integer *, doublereal *, integer *, integer *, ftnlen);
  static logical notran;
  static char transt[1];
  static T lstres;

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  afb_dim1 = *ldafb;
  afb_offset = 1 + afb_dim1;
  afb -= afb_offset;
  --ipiv;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;
  x_dim1 = *ldx;
  x_offset = 1 + x_dim1;
  x -= x_offset;
  --ferr;
  --berr;
  --work;
  --iwork;

  /* Function Body */
  *info = 0;
  notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
  if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
    trans, "C", (ftnlen)1, (ftnlen)1)) {
    *info = -1;
  } else if (*n < 0) {
    *info = -2;
  } else if (*kl < 0) {
    *info = -3;
  } else if (*ku < 0) {
    *info = -4;
  } else if (*nrhs < 0) {
    *info = -5;
  } else if (*ldab < *kl + *ku + 1) {
    *info = -7;
  } else if (*ldafb < (*kl << 1) + *ku + 1) {
    *info = -9;
  } else if (*ldb < max(1,*n)) {
    *info = -12;
  } else if (*ldx < max(1,*n)) {
    *info = -14;
  }
  if (*info != 0) {
    i__1 = -(*info);
    xerbla_("DGBRFS", &i__1, (ftnlen)6);
    return 0;
  }

  /*     Quick return if possible */

  if (*n == 0 || *nrhs == 0) {
    i__1 = *nrhs;
    for (j = 1; j <= i__1; ++j) {
      ferr[j] = 0.;
      berr[j] = 0.;
      /* L10: */
    }
    return 0;
  }

  if (notran) {
    *(unsigned char *)transt = 'T';
  } else {
    *(unsigned char *)transt = 'N';
  }

  /*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

  /* Computing MIN */
  i__1 = *kl + *ku + 2, i__2 = *n + 1;
  nz = min(i__1,i__2);
  eps = dlamch_("Epsilon", (ftnlen)7);
  safmin = dlamch_("Safe minimum", (ftnlen)12);
  safe1 = nz * safmin;
  safe2 = safe1 / eps;

  /*     Do for each right hand side */

  i__1 = *nrhs;
  for (j = 1; j <= i__1; ++j) {

    count = 1;
    lstres = 3.;
  L20:

    /*        Loop until stopping criterion is satisfied. */

    /*        Compute residual R = B - op(A) * X, */
    /*        where op(A) = A, A**T, or A**H, depending on TRANS. */

    dcopy_(n, &b[j * b_dim1 + 1], &c__1, &work[*n + 1], &c__1);
    dgbmv_(trans, n, n, kl, ku, &c_b15, &ab[ab_offset], ldab, &x[j * 
      x_dim1 + 1], &c__1, &c_b17, &work[*n + 1], &c__1, (ftnlen)1);

    /*        Compute componentwise relative backward error from formula */

    /*        max(i) ( abs(R(i)) / ( abs(op(A))*abs(X) + abs(B) )(i) ) */

    /*        where abs(Z) is the componentwise absolute value of the matrix */
    /*        or vector Z.  If the i-th component of the denominator is less */
    /*        than SAFE2, then SAFE1 is added to the i-th components of the */
    /*        numerator and denominator before dividing. */

    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      work[i__] = (d__1 = b[i__ + j * b_dim1], abs(d__1));
      /* L30: */
    }

    /*        Compute abs(op(A))*abs(X) + abs(B). */

    if (notran) {
      i__2 = *n;
      for (k = 1; k <= i__2; ++k) {
        kk = *ku + 1 - k;
        xk = (d__1 = x[k + j * x_dim1], abs(d__1));
        /* Computing MAX */
        i__3 = 1, i__4 = k - *ku;
        /* Computing MIN */
        i__6 = *n, i__7 = k + *kl;
        i__5 = min(i__6,i__7);
        for (i__ = max(i__3,i__4); i__ <= i__5; ++i__) {
          work[i__] += (d__1 = ab[kk + i__ + k * ab_dim1], abs(d__1)
            ) * xk;
          /* L40: */
        }
        /* L50: */
      }
    } else {
      i__2 = *n;
      for (k = 1; k <= i__2; ++k) {
        s = 0.;
        kk = *ku + 1 - k;
        /* Computing MAX */
        i__5 = 1, i__3 = k - *ku;
        /* Computing MIN */
        i__6 = *n, i__7 = k + *kl;
        i__4 = min(i__6,i__7);
        for (i__ = max(i__5,i__3); i__ <= i__4; ++i__) {
          s += (d__1 = ab[kk + i__ + k * ab_dim1], abs(d__1)) * (
            d__2 = x[i__ + j * x_dim1], abs(d__2));
          /* L60: */
        }
        work[k] += s;
        /* L70: */
      }
    }
    s = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      if (work[i__] > safe2) {
        /* Computing MAX */
        d__2 = s, d__3 = (d__1 = work[*n + i__], abs(d__1)) / work[
          i__];
        s = max(d__2,d__3);
      } else {
        /* Computing MAX */
        d__2 = s, d__3 = ((d__1 = work[*n + i__], abs(d__1)) + safe1) 
          / (work[i__] + safe1);
        s = max(d__2,d__3);
      }
      /* L80: */
    }
    berr[j] = s;

    /*        Test stopping criterion. Continue iterating if */
    /*           1) The residual BERR(J) is larger than machine epsilon, and */
    /*           2) BERR(J) decreased by at least a factor of 2 during the */
    /*              last iteration, and */
    /*           3) At most ITMAX iterations tried. */

    if (berr[j] > eps && berr[j] * 2. <= lstres && count <= 5) {

      /*           Update solution and try again. */

      dgbtrs_(trans, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &ipiv[1]
        , &work[*n + 1], n, info, (ftnlen)1);
      daxpy_(n, &c_b17, &work[*n + 1], &c__1, &x[j * x_dim1 + 1], &c__1)
        ;
      lstres = berr[j];
      ++count;
      goto L20;
    }

    /*        Bound error from formula */

    /*        norm(X - XTRUE) / norm(X) .le. FERR = */
    /*        norm( abs(inv(op(A)))* */
    /*           ( abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) / norm(X) */

    /*        where */
    /*          norm(Z) is the magnitude of the largest component of Z */
    /*          inv(op(A)) is the inverse of op(A) */
    /*          abs(Z) is the componentwise absolute value of the matrix or */
    /*             vector Z */
    /*          NZ is the maximum number of nonzeros in any row of A, plus 1 */
    /*          EPS is machine epsilon */

    /*        The i-th component of abs(R)+NZ*EPS*(abs(op(A))*abs(X)+abs(B)) */
    /*        is incremented by SAFE1 if the i-th component of */
    /*        abs(op(A))*abs(X) + abs(B) is less than SAFE2. */

    /*        Use DLACN2 to estimate the infinity-norm of the matrix */
    /*           inv(op(A)) * diag(W), */
    /*        where W = abs(R) + NZ*EPS*( abs(op(A))*abs(X)+abs(B) ))) */

    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      if (work[i__] > safe2) {
        work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
          work[i__];
      } else {
        work[i__] = (d__1 = work[*n + i__], abs(d__1)) + nz * eps * 
          work[i__] + safe1;
      }
      /* L90: */
    }

    kase = 0;
  L100:
    dlacn2_(n, &work[(*n << 1) + 1], &work[*n + 1], &iwork[1], &ferr[j], &
      kase, isave);
    if (kase != 0) {
      if (kase == 1) {

        /*              Multiply by diag(W)*inv(op(A)**T). */

        dgbtrs_(transt, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &
          ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          work[*n + i__] *= work[i__];
          /* L110: */
        }
      } else {

        /*              Multiply by inv(op(A))*diag(W). */

        i__2 = *n;
        for (i__ = 1; i__ <= i__2; ++i__) {
          work[*n + i__] *= work[i__];
          /* L120: */
        }
        dgbtrs_(trans, n, kl, ku, &c__1, &afb[afb_offset], ldafb, &
          ipiv[1], &work[*n + 1], n, info, (ftnlen)1);
      }
      goto L100;
    }

    /*        Normalize error. */

    lstres = 0.;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
      /* Computing MAX */
      d__2 = lstres, d__3 = (d__1 = x[i__ + j * x_dim1], abs(d__1));
      lstres = max(d__2,d__3);
      /* L130: */
    }
    if (lstres != 0.) {
      ferr[j] /= lstres;
    }

    /* L140: */
  }

  return 0;

  /*     End of DGBRFS */

} /* dgbrfs_ */

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

// translated from dlangb - LAPACK auxiliary routine (version 3.7.0) */
template<class T>
doublereal langb(char *norm, integer *n, integer *kl, integer *ku, T *ab, integer *ldab, T *work, 
  ftnlen norm_len)
{
  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;
  T ret_val, d__1;

  /* Builtin functions */
  double sqrt(T);

  /* Local variables */
  static integer i__, j, k, l;
  static T sum, temp, scale;
  extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static T value;
  extern logical disnan_(doublereal *);
  extern /* Subroutine */ int dlassq_(integer *, T *, integer *, T *, T *);

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --work;

  /* Function Body */
  if (*n == 0) {
    value = 0.;
  } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

    /*        Find max(abs(A(i,j))). */

    value = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      /* Computing MAX */
      i__2 = *ku + 2 - j;
      /* Computing MIN */
      i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
      i__3 = min(i__4,i__5);
      for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
        temp = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
        if (value < temp || disnan_(&temp)) {
          value = temp;
        }
        /* L10: */
      }
      /* L20: */
    }
  } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
    norm == '1') {

    /*        Find norm1(A). */

    value = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      sum = 0.;
      /* Computing MAX */
      i__3 = *ku + 2 - j;
      /* Computing MIN */
      i__4 = *n + *ku + 1 - j, i__5 = *kl + *ku + 1;
      i__2 = min(i__4,i__5);
      for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
        sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
        /* L30: */
      }
      if (value < sum || disnan_(&sum)) {
        value = sum;
      }
      /* L40: */
    }
  } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

    /*        Find normI(A). */

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      work[i__] = 0.;
      /* L50: */
    }
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      k = *ku + 1 - j;
      /* Computing MAX */
      i__2 = 1, i__3 = j - *ku;
      /* Computing MIN */
      i__5 = *n, i__6 = j + *kl;
      i__4 = min(i__5,i__6);
      for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
        work[i__] += (d__1 = ab[k + i__ + j * ab_dim1], abs(d__1));
        /* L60: */
      }
      /* L70: */
    }
    value = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      temp = work[i__];
      if (value < temp || disnan_(&temp)) {
        value = temp;
      }
      /* L80: */
    }
  } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
    ftnlen)1, (ftnlen)1)) {

    /*        Find normF(A). */

    scale = 0.;
    sum = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      /* Computing MAX */
      i__4 = 1, i__2 = j - *ku;
      l = max(i__4,i__2);
      k = *ku + 1 - j + l;
      /* Computing MIN */
      i__2 = *n, i__3 = j + *kl;
      i__4 = min(i__2,i__3) - l + 1;
      dlassq_(&i__4, &ab[k + j * ab_dim1], &c__1, &scale, &sum);
      /* L90: */
    }
    value = scale * sqrt(sum);
  }

  ret_val = value;
  return ret_val;

  /*     End of DLANGB */

} /* dlangb_ */

//-------------------------------------------------------------------------------------------------

// translated from dlacpy - LAPACK auxiliary routine (version 3.7.0) 
template<class T>
int lacpy(char *uplo, integer *m, integer *n, T *a, integer *lda, T *b, integer *ldb, 
  ftnlen uplo_len)
{
  /* System generated locals */
  integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

  /* Local variables */
  static integer i__, j;
  extern logical lsame_(char *, char *, ftnlen, ftnlen);

  /* Parameter adjustments */
  a_dim1 = *lda;
  a_offset = 1 + a_dim1;
  a -= a_offset;
  b_dim1 = *ldb;
  b_offset = 1 + b_dim1;
  b -= b_offset;

  /* Function Body */
  if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      i__2 = min(j,*m);
      for (i__ = 1; i__ <= i__2; ++i__) {
        b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
        /* L10: */
      }
      /* L20: */
    }
  } else if (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      i__2 = *m;
      for (i__ = j; i__ <= i__2; ++i__) {
        b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
        /* L30: */
      }
      /* L40: */
    }
  } else {
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      i__2 = *m;
      for (i__ = 1; i__ <= i__2; ++i__) {
        b[i__ + j * b_dim1] = a[i__ + j * a_dim1];
        /* L50: */
      }
      /* L60: */
    }
  }
  return 0;

} /* dlacpy_ */

//-------------------------------------------------------------------------------------------------

// from dlantb - LAPACK auxiliary routine (version 3.7.0)
template<class T>
T lantb(char *norm, char *uplo, char *diag, integer *n, integer *k,
  T *ab, integer *ldab, doublereal *work, ftnlen norm_len, 
  ftnlen uplo_len, ftnlen diag_len)
{
  /* Table of constant values */
  static integer c__1 = 1;

  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5;
  T ret_val, d__1;

  /* Builtin functions */
  double sqrt(T);

  /* Local variables */
  static integer i__, j, l;
  static T sum, scale;
  static logical udiag;
  extern logical lsame_(char *, char *, ftnlen, ftnlen);
  static T value;
  extern logical disnan_(T *);
  extern /* Subroutine */ int dlassq_(integer *, T *, integer *, T *, T *);

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --work;

  /* Function Body */
  if (*n == 0) {
    value = 0.;
  } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

    /*        Find max(abs(A(i,j))). */

    if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
      value = 1.;
      if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MAX */
          i__2 = *k + 2 - j;
          i__3 = *k;
          for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
            sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            if (value < sum || disnan_(&sum)) {
              value = sum;
            }
            /* L10: */
          }
          /* L20: */
        }
      } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MIN */
          i__2 = *n + 1 - j, i__4 = *k + 1;
          i__3 = min(i__2,i__4);
          for (i__ = 2; i__ <= i__3; ++i__) {
            sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            if (value < sum || disnan_(&sum)) {
              value = sum;
            }
            /* L30: */
          }
          /* L40: */
        }
      }
    } else {
      value = 0.;
      if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MAX */
          i__3 = *k + 2 - j;
          i__2 = *k + 1;
          for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
            sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            if (value < sum || disnan_(&sum)) {
              value = sum;
            }
            /* L50: */
          }
          /* L60: */
        }
      } else {
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MIN */
          i__3 = *n + 1 - j, i__4 = *k + 1;
          i__2 = min(i__3,i__4);
          for (i__ = 1; i__ <= i__2; ++i__) {
            sum = (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            if (value < sum || disnan_(&sum)) {
              value = sum;
            }
            /* L70: */
          }
          /* L80: */
        }
      }
    }
  } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
    norm == '1') {

    /*        Find norm1(A). */

    value = 0.;
    udiag = lsame_(diag, "U", (ftnlen)1, (ftnlen)1);
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        if (udiag) {
          sum = 1.;
          /* Computing MAX */
          i__2 = *k + 2 - j;
          i__3 = *k;
          for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
            sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            /* L90: */
          }
        } else {
          sum = 0.;
          /* Computing MAX */
          i__3 = *k + 2 - j;
          i__2 = *k + 1;
          for (i__ = max(i__3,1); i__ <= i__2; ++i__) {
            sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            /* L100: */
          }
        }
        if (value < sum || disnan_(&sum)) {
          value = sum;
        }
        /* L110: */
      }
    } else {
      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        if (udiag) {
          sum = 1.;
          /* Computing MIN */
          i__3 = *n + 1 - j, i__4 = *k + 1;
          i__2 = min(i__3,i__4);
          for (i__ = 2; i__ <= i__2; ++i__) {
            sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            /* L120: */
          }
        } else {
          sum = 0.;
          /* Computing MIN */
          i__3 = *n + 1 - j, i__4 = *k + 1;
          i__2 = min(i__3,i__4);
          for (i__ = 1; i__ <= i__2; ++i__) {
            sum += (d__1 = ab[i__ + j * ab_dim1], abs(d__1));
            /* L130: */
          }
        }
        if (value < sum || disnan_(&sum)) {
          value = sum;
        }
        /* L140: */
      }
    }
  } else if (lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

    /*        Find normI(A). */

    value = 0.;
    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
      if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] = 1.;
          /* L150: */
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          l = *k + 1 - j;
          /* Computing MAX */
          i__2 = 1, i__3 = j - *k;
          i__4 = j - 1;
          for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
            work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
              d__1));
            /* L160: */
          }
          /* L170: */
        }
      } else {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] = 0.;
          /* L180: */
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          l = *k + 1 - j;
          /* Computing MAX */
          i__4 = 1, i__2 = j - *k;
          i__3 = j;
          for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
            work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
              d__1));
            /* L190: */
          }
          /* L200: */
        }
      }
    } else {
      if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] = 1.;
          /* L210: */
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          l = 1 - j;
          /* Computing MIN */
          i__4 = *n, i__2 = j + *k;
          i__3 = min(i__4,i__2);
          for (i__ = j + 1; i__ <= i__3; ++i__) {
            work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
              d__1));
            /* L220: */
          }
          /* L230: */
        }
      } else {
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          work[i__] = 0.;
          /* L240: */
        }
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          l = 1 - j;
          /* Computing MIN */
          i__4 = *n, i__2 = j + *k;
          i__3 = min(i__4,i__2);
          for (i__ = j; i__ <= i__3; ++i__) {
            work[i__] += (d__1 = ab[l + i__ + j * ab_dim1], abs(
              d__1));
            /* L250: */
          }
          /* L260: */
        }
      }
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
      sum = work[i__];
      if (value < sum || disnan_(&sum)) {
        value = sum;
      }
      /* L270: */
    }
  } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
    ftnlen)1, (ftnlen)1)) {

    /*        Find normF(A). */

    if (lsame_(uplo, "U", (ftnlen)1, (ftnlen)1)) {
      if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
        scale = 1.;
        sum = (doublereal) (*n);
        if (*k > 0) {
          i__1 = *n;
          for (j = 2; j <= i__1; ++j) {
            /* Computing MIN */
            i__4 = j - 1;
            i__3 = min(i__4,*k);
            /* Computing MAX */
            i__2 = *k + 2 - j;
            dlassq_(&i__3, &ab[max(i__2,1) + j * ab_dim1], &c__1, 
              &scale, &sum);
            /* L280: */
          }
        }
      } else {
        scale = 0.;
        sum = 1.;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MIN */
          i__4 = j, i__2 = *k + 1;
          i__3 = min(i__4,i__2);
          /* Computing MAX */
          i__5 = *k + 2 - j;
          dlassq_(&i__3, &ab[max(i__5,1) + j * ab_dim1], &c__1, &
            scale, &sum);
          /* L290: */
        }
      }
    } else {
      if (lsame_(diag, "U", (ftnlen)1, (ftnlen)1)) {
        scale = 1.;
        sum = (doublereal) (*n);
        if (*k > 0) {
          i__1 = *n - 1;
          for (j = 1; j <= i__1; ++j) {
            /* Computing MIN */
            i__4 = *n - j;
            i__3 = min(i__4,*k);
            dlassq_(&i__3, &ab[j * ab_dim1 + 2], &c__1, &scale, &
              sum);
            /* L300: */
          }
        }
      } else {
        scale = 0.;
        sum = 1.;
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          /* Computing MIN */
          i__4 = *n - j + 1, i__2 = *k + 1;
          i__3 = min(i__4,i__2);
          dlassq_(&i__3, &ab[j * ab_dim1 + 1], &c__1, &scale, &sum);
          /* L310: */
        }
      }
    }
    value = scale * sqrt(sum);
  }

  ret_val = value;
  return ret_val;

  /*     End of DLANTB */

} /* dlantb_ */

//-------------------------------------------------------------------------------------------------

// from dlaqgb - LAPACK auxiliary routine (version 3.7.0) 
template<class T>
int laqgb(integer *m, integer *n, integer *kl, integer *ku,
  T *ab, integer *ldab, T *r__, T *c__, T *rowcnd, T *colcnd, T *amax, char *equed, 
  ftnlen equed_len)
{
  /* System generated locals */
  integer ab_dim1, ab_offset, i__1, i__2, i__3, i__4, i__5, i__6;

  /* Local variables */
  static integer i__, j;
  static T cj, large, small;
  extern T dlamch_(char *, ftnlen);

  /* Parameter adjustments */
  ab_dim1 = *ldab;
  ab_offset = 1 + ab_dim1;
  ab -= ab_offset;
  --r__;
  --c__;

  /* Function Body */
  if (*m <= 0 || *n <= 0) {
    *(unsigned char *)equed = 'N';
    return 0;
  }

  /*     Initialize LARGE and SMALL. */

  small = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
    ftnlen)9);
  large = 1. / small;

  if (*rowcnd >= .1 && *amax >= small && *amax <= large) {

    /*        No row scaling */

    if (*colcnd >= .1) {

      /*           No column scaling */

      *(unsigned char *)equed = 'N';
    } else {

      /*           Column scaling */

      i__1 = *n;
      for (j = 1; j <= i__1; ++j) {
        cj = c__[j];
        /* Computing MAX */
        i__2 = 1, i__3 = j - *ku;
        /* Computing MIN */
        i__5 = *m, i__6 = j + *kl;
        i__4 = min(i__5,i__6);
        for (i__ = max(i__2,i__3); i__ <= i__4; ++i__) {
          ab[*ku + 1 + i__ - j + j * ab_dim1] = cj * ab[*ku + 1 + 
            i__ - j + j * ab_dim1];
          /* L10: */
        }
        /* L20: */
      }
      *(unsigned char *)equed = 'C';
    }
  } else if (*colcnd >= .1) {

    /*        Row scaling, no column scaling */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      /* Computing MAX */
      i__4 = 1, i__2 = j - *ku;
      /* Computing MIN */
      i__5 = *m, i__6 = j + *kl;
      i__3 = min(i__5,i__6);
      for (i__ = max(i__4,i__2); i__ <= i__3; ++i__) {
        ab[*ku + 1 + i__ - j + j * ab_dim1] = r__[i__] * ab[*ku + 1 + 
          i__ - j + j * ab_dim1];
        /* L30: */
      }
      /* L40: */
    }
    *(unsigned char *)equed = 'R';
  } else {

    /*        Row and column scaling */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
      cj = c__[j];
      /* Computing MAX */
      i__3 = 1, i__4 = j - *ku;
      /* Computing MIN */
      i__5 = *m, i__6 = j + *kl;
      i__2 = min(i__5,i__6);
      for (i__ = max(i__3,i__4); i__ <= i__2; ++i__) {
        ab[*ku + 1 + i__ - j + j * ab_dim1] = cj * r__[i__] * ab[*ku 
          + 1 + i__ - j + j * ab_dim1];
        /* L50: */
      }
      /* L60: */
    }
    *(unsigned char *)equed = 'B';
  }

  return 0;

  /*     End of DLAQGB */

} /* dlaqgb_ */

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

template int gbsv(long int *n, long int *kl, long int *ku, long int *nrhs, double *ab, 
  long int *ldab, long int *ipiv, double *b, long int *ldb, long int *info);


}