#pragma once

/** For some of the functions in LibF2C, Blas and LaPack, it didn't make sense to templatize them.
These are all put into this file because it may be problematic in certain contexts to have mixed 
templated and typed code in a single implementation file - if you include the same implementation
file in two different unity-build compilation units, you will get "multiple definition" linker 
errors. You may have to include the same implementation file into different compilation unit, if
these compilation units instantiate different templates or the same template for different 
datatypes. I know - it's a mess :-( */

namespace LaPackCPP {

//=================================================================================================
// LibF2C

integer i_len(char *s, ftnlen n)
{
  return(n);
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


//=================================================================================================
// Blas


//=================================================================================================
// XBlas

void BLAS_error(const char *rname, int iflag, int ival, char *form, ...)
{
#if !defined(CONFIG_USE_XERBLA)
{
  va_list argptr;
  va_start(argptr, form);
  fprintf(stderr, "Error #%d from routine %s:\n", iflag, rname);
  if(form)
    vfprintf(stderr, form, argptr);
  else if(iflag < 0)
    fprintf(stderr,
      "  Parameter number %d to routine %s had the illegal value %d\n",
      -iflag, rname, ival);
  else
    fprintf(stderr, "  Unknown error code %d from routine %s\n",
      iflag, rname);
  exit(iflag);
}
#else
{
  int ln, argno;
  ln = strlen(rname);
  argno = -iflag;
  xerbla_array(rname, &ln, &argno);
}
#endif
}

//=================================================================================================
// LaPack







}