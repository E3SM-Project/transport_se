/* $Id: timer.c,v 2.2 2005/07/19 16:34:54 mataylo Exp $ */

#if defined(Darwin)
#include <sys/time.h>
double timer_clock()
{
  double t;
  struct timeval buffer;
  struct timezone dummy;
  gettimeofday (&buffer, &dummy);
  t = (double)(buffer.tv_sec*1000000 + buffer.tv_usec);
  return (t);
}
#endif

#if defined(__hpux)
#include <time.h>
#include <stdio.h>
#include <math.h>
double timer_clock()
{
  double t;
  struct timeval buffer;
  struct timezone dummy;
  gettimeofday (&buffer, &dummy);
  t = (double)(buffer.tv_sec*1000000 + buffer.tv_usec);
  return (t);
}

#elif defined (__sun)
#include <sys/time.h>
static int first_call = 1;
double timer_clock()
{
  hrtime_t nsec;
  /*  if (first_call) 
      {
      first_call = 0;
      init_ecache_();
      } 
  */
  nsec = gethrtime();
  return((double)nsec*1.0e-03);
}

#elif defined(Linux) 
#include <sys/time.h>
double timer_clock()
{
  double t;
  struct timeval buffer;
  struct timezone dummy;
  gettimeofday (&buffer, &dummy);
  t = (double)(buffer.tv_sec*1000000 + buffer.tv_usec);
  return (t);
}

#elif defined(TFLOPS) 
#include <sys/time.h>
double timer_clock()
{
#include <nx.h>
  double t;
  t=1e6*dclock();
  return (t);
}

#elif defined(__alpha) 
#include <sys/time.h>
double timer_clock()
{
  double t;
  struct timeval buffer;
  struct timezone dummy;
  gettimeofday (&buffer, &dummy);
  t = (double)(buffer.tv_sec*1000000+buffer.tv_usec);
  return (t);
}

#elif defined(__sgi) 
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
double timer_clock()
{
  struct rusage 
  rusage1, rusage2;
  getrusage(RUSAGE_SELF, &rusage2);
  return( (double) (rusage2.ru_utime.tv_sec*1000000 +
		    rusage2.ru_utime.tv_usec  +
		    rusage2.ru_stime.tv_sec*1000000 +
		    rusage2.ru_stime.tv_usec));
}

#elif defined(AIX) 
double rtc(void); /* note: must link with -lxlf90 to get rtc() */
double timer_clock()
{
  double t;
  t = rtc()*1.0e6;
  return (t);
}
#endif

/* Fortran callable wrappers */

double timer_clock__(){
double timer_clock();
return(timer_clock());
}
double timer_clock_(){
double timer_clock();
return(timer_clock());
}

double TIMER_CLOCK(){
double timer_clock();
return(timer_clock());
}
