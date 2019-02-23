#include "f2c.h"
#ifdef __cplusplus
extern "C" {
#endif

#ifdef KR_headers
double r_dim(a,b) real *a, *b;
#else
double r_dim(f2c_real *a, f2c_real *b)
#endif
{
return( *a > *b ? *a - *b : 0);
}
#ifdef __cplusplus
}
#endif
