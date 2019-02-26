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

// translated from lsame, Reference BLAS level1 routine (version 3.1)
logical lsame(char *ca, char *cb, ftnlen ca_len, ftnlen cb_len)
{
  // System generated locals
  logical ret_val;

  // Local variables 
  static integer inta, intb, zcode;

  // Test if the characters are equal
  ret_val = *(unsigned char *)ca == *(unsigned char *)cb;
  if(ret_val) {
    return ret_val;
  }

  // Now test for equivalence if both characters are alphabetic.
  zcode = 'Z';

  // Use 'Z' rather than 'A' so that ASCII can be detected on Prime 
  // machines, on which ICHAR returns a value with bit 8 set. 
  // ICHAR('A') on Prime machines returns 193 which is the same as 
  // ICHAR('A') on an EBCDIC machine.

  inta = *(unsigned char *)ca;
  intb = *(unsigned char *)cb;

  if(zcode == 90 || zcode == 122) {

    // ASCII is assumed - ZCODE is the ASCII code of either lower or upper case 'Z'. 
    if(inta >= 97 && inta <= 122) {
      inta += -32;
    }
    if(intb >= 97 && intb <= 122) {
      intb += -32;
    }

  }
  else if(zcode == 233 || zcode == 169) {

    // EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or 
    // upper case 'Z'.
    if(inta >= 129 && inta <= 137 || inta >= 145 && inta <= 153 || inta
      >= 162 && inta <= 169) {
      inta += 64;
    }
    if(intb >= 129 && intb <= 137 || intb >= 145 && intb <= 153 || intb
      >= 162 && intb <= 169) {
      intb += 64;
    }

  }
  else if(zcode == 218 || zcode == 250) {

    // ASCII is assumed, on Prime machines - ZCODE is the ASCII code
    // plus 128 of either lower or upper case 'Z'. 

    if(inta >= 225 && inta <= 250) {
      inta += -32;
    }
    if(intb >= 225 && intb <= 250) {
      intb += -32;
    }
  }
  ret_val = inta == intb;

  // End of LSAME

  return ret_val;
}


int xerbla(char *srname, integer *info, ftnlen srname_len)
{
  // Code of the function has been commented out by Robin Schmidt - at the moment, xerbla is just
  // a dummy function - i should probably set a debug-breakpoint here and re-implement it 
  // completely - it deals with some weird I/O functions from the f2c library which seems to be 
  // useless clutter...

  /*
  // Table of constant values
  static integer c__1 = 1;

  // Format strings 
  static char fmt_9999[] = "(\002 ** On entry to \002,a,\002 parameter num"
    "ber \002,i2,\002 had \002,\002an illegal value\002)";

  // Builtin functions 
  integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
  int s_stop(char *, ftnlen);

  // Local variables
  extern integer len_trim__(char *, ftnlen);
  // Note by Robin Schmidt:
  // this was declared as an intrinsic function in the original xerbla.f file, but it stumped the 
  // f2c translator, giving an error about an unknown intrinsic function, so i changed the 
  // "INTRINSIC" keyword to "EXTERNAL" - that allowed the translation, but i'm not sure, if it 
  // gives the correct behavior 

  // Fortran I/O blocks
  static cilist io___1 = { 0, 6, 0, fmt_9999, 0 };

  s_wsfe(&io___1);
  do_fio(&c__1, srname, len_trim__(srname, srname_len));
  do_fio(&c__1, (char *)&(*info), (ftnlen)sizeof(integer));
  e_wsfe();

  s_stop("", (ftnlen)0);

  // End of XERBLA
  */

  return 0;
} // xerbla


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