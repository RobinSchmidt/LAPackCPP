#include "f2c.h"

#ifdef KR_headers
double atan2();
double r_atn2(x,y) real *x, *y;
#else
#undef abs
#include "math.h"
#ifdef __cplusplus
extern "C" {
#endif
double r_atn2(f2c_real *x, f2c_real *y)
#endif
{
return( atan2(*x,*y) );
}
#ifdef __cplusplus
}
#endif
