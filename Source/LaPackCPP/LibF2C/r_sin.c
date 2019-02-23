#include "f2c.h"

#ifdef KR_headers
double sin();
double r_sin(x) real *x;
#else
#undef abs
#include "math.h"
#ifdef __cplusplus
extern "C" {
#endif
double r_sin(f2c_real *x)
#endif
{
return( sin(*x) );
}
#ifdef __cplusplus
}
#endif
