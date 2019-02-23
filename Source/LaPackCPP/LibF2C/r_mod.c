#include "f2c.h"

#ifdef KR_headers
#ifdef IEEE_drem
double drem();
#else
double floor();
#endif
double r_mod(x,y) real *x, *y;
#else
#ifdef IEEE_drem
double drem(double, double);
#else
#undef abs
#include "math.h"
#ifdef __cplusplus
extern "C" {
#endif
#endif
double r_mod(f2c_real *x, f2c_real *y)
#endif
{
#ifdef IEEE_drem
	double xa, ya, z;
	if ((ya = *y) < 0.)
		ya = -ya;
	z = drem(xa = *x, ya);
	if (xa > 0) {
		if (z < 0)
			z += ya;
		}
	else if (z > 0)
		z -= ya;
	return z;
#else
	double quotient;
	if( (quotient = (double)*x / *y) >= 0)
		quotient = floor(quotient);
	else
		quotient = -floor(-quotient);
	return(*x - (*y) * quotient );
#endif
}
#ifdef __cplusplus
}
#endif
