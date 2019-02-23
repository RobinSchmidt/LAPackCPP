#include "f2c.h"

#ifdef KR_headers
double log();
double r_log(x) real *x;
#else
#undef abs
#include "math.h"
#ifdef __cplusplus
extern "C" {
#endif
double r_log(f2c_real *x)
#endif
{
return( log(*x) );
}
#ifdef __cplusplus
}
#endif
