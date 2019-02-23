#include <cmath>
#include "LibF2C.hpp"


#define log10e 0.43429448190325182765

namespace LibF2C {

// todo: templatize all these functions

double d_lg10(doublereal *x)
{
  return( log10e * log(*x) );
}

double d_sign(doublereal *a, doublereal *b)
{
  double x;
  x = (*a >= 0 ? *a : - *a);
  return( *b >= 0 ? x : -x);
}

integer i_dnnt(doublereal *x)
{
  return (integer)(*x >= 0. ? floor(*x + .5) : -floor(.5 - *x));
}

integer i_len(char *s, ftnlen n)
{
  return(n);
}

integer i_nint(f2c_real *x)
{
  return (integer)(*x >= 0 ? floor(*x + .5) : -floor(.5 - *x));
}

double pow_di(doublereal *ap, integer *bp)
{
  double pow, x;
  integer n;
  unsigned long u;

  pow = 1;
  x = *ap;
  n = *bp;

  if(n != 0)
  {
    if(n < 0)
    {
      n = -n;
      x = 1/x;
    }
    for(u = n; ; )
    {
      if(u & 01)
        pow *= x;
      if(u >>= 1)
        x *= x;
      else
        break;
    }
  }
  return(pow);
}

integer s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
{
  register unsigned char *a, *aend, *b, *bend;
  a = (unsigned char *)a0;
  b = (unsigned char *)b0;
  aend = a + la;
  bend = b + lb;

  if(la <= lb)
  {
    while(a < aend)
      if(*a != *b)
        return( *a - *b );
      else
      { ++a; ++b; }

    while(b < bend)
      if(*b != ' ')
        return( ' ' - *b );
      else	++b;
  }

  else
  {
    while(b < bend)
      if(*a == *b)
      { ++a; ++b; }
      else
        return( *a - *b );
    while(a < aend)
      if(*a != ' ')
        return(*a - ' ');
      else	++a;
  }
  return(0);
}

void s_copy(register char *a, register char *b, ftnlen la, ftnlen lb)
{
  register char *aend, *bend;

  aend = a + la;

  if(la <= lb)
#ifndef NO_OVERWRITE
    if(a <= b || a >= b + la)
#endif
      while(a < aend)
        *a++ = *b++;
#ifndef NO_OVERWRITE
    else
      for(b += la; a < aend; )
        *--aend = *--b;
#endif

  else {
    bend = b + lb;
#ifndef NO_OVERWRITE
    if(a <= b || a >= bend)
#endif
      while(b < bend)
        *a++ = *b++;
#ifndef NO_OVERWRITE
    else {
      a += lb;
      while(b < bend)
        *--a = *--bend;
      a += lb;
    }
#endif
    while(a < aend)
      *a++ = ' ';
  }
}




}