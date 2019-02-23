#include "time.h"

#ifdef MSDOS
#undef USE_CLOCK
#define USE_CLOCK
#endif

#ifndef REAL
#define REAL double
#endif

#ifndef USE_CLOCK
#define _INCLUDE_POSIX_SOURCE	/* for HP-UX */
#define _INCLUDE_XOPEN_SOURCE	/* for HP-UX */
#include "sys/types.h"

//#include "sys/times.h" // Robin Schmidt: this doesn't compile in MSVC, infos here
// http://pubs.opengroup.org/onlinepubs/7908799/xsh/systimes.h.html
// so i just added the definition of the tms struct directly here from here:
// https://code.woboq.org/gcc/include/sys/times.h.html
// instead of including sys/times.h:
/* Structure describing CPU time used by a process and its children.  */
struct tms
{
  clock_t tms_utime;                /* User CPU time.  */
  clock_t tms_stime;                /* System CPU time.  */
  clock_t tms_cutime;               /* User CPU time of dead children.  */
  clock_t tms_cstime;               /* System CPU time of dead children.  */
};


#ifdef __cplusplus
extern "C" {
#endif
#endif

#undef Hz
#ifdef CLK_TCK
#define Hz CLK_TCK
#else
#ifdef HZ
#define Hz HZ
#else
#define Hz 60
#endif
#endif

 REAL
#ifdef KR_headers
etime_(tarray) float *tarray;
#else
etime_(float *tarray)
#endif
{
#ifdef USE_CLOCK
#ifndef CLOCKS_PER_SECOND
#define CLOCKS_PER_SECOND Hz
#endif
	double t = clock();
	tarray[1] = 0;
	return tarray[0] = t / CLOCKS_PER_SECOND;
#else
	struct tms t;

	times(&t);
	return	  (tarray[0] = (double)t.tms_utime/Hz)
		+ (tarray[1] = (double)t.tms_stime/Hz);
#endif
	}
#ifdef __cplusplus
}
#endif
