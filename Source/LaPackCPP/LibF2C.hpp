#pragma once

namespace LibF2C {

typedef long int integer;
//typedef unsigned long int uinteger; // added by Robin Schmidt (guesswork!!)
//typedef char *address;
//typedef short int shortint;
typedef float f2c_real;
typedef double doublereal;
//typedef struct { f2c_real r, i; } f2c_complex;
//typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
//typedef short int shortlogical;
typedef long ftnlen;

#define TRUE_ (1)
#define FALSE_ (0)
#define VOID void

// try to get rid of them:
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))


extern double d_lg10(double *);
extern double d_sign(double *, double *);
extern integer i_dnnt(double *);
extern integer i_len(char *, ftnlen);
extern integer i_nint(float *);
extern double pow_di(double *, integer *);
extern integer s_cmp(char *, char *, ftnlen, ftnlen);
extern void s_copy(char *, char *, ftnlen, ftnlen);


}