#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifndef REAL
#define REAL double
#endif

#ifdef KR_headers
double erfc();
REAL erfc_(x) real *x;
#else
extern double erfc(double);
REAL erfc_(f2c_real *x)
#endif
{
return( erfc((double)*x) );
}
#ifdef __cplusplus
}
#endif
