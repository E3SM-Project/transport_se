/*
   =================================================
   Compute the number of bits of precision in a double
   floating point type 
   =================================================
*/

#if defined(Linux) || defined(SunOS) || defined(OSF1) || defined(TFLOPS)
int precision_(){
#endif
#if defined(AIX)
int precision(){
#endif
#if defined(Darwin)
int precision_(){
#endif

double a=1.0;
double a2=1.0;
double one=1.0;
int iprec=0;
while(a+a2 != one){
 a2=a2/2;
 iprec++;
}
iprec--;

return(iprec);

}

/*
   ===========================================================
   Compute the number of bits of precision in a long double
   floating point type 
   ===========================================================
*/

#if defined(Linux) || defined(SunOS) || defined(OSF1) || defined(TFLOPS)
int dprecision_(){
#endif
#if defined(AIX)
int dprecision(){
#endif
#if defined(Darwin)
int dprecision_(){
#endif

long double a=1.0;
long double a2=1.0;
long double one=1.0;
long double two=2.0;
int iprec=0;
while(a+a2 != one){
 a2=a2/two;
 iprec++;
}
iprec--;

return(iprec);

}
