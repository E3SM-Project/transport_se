/* Computes the number of bits needed to
   represent a number */

int nbits(int n){
  int nh=1;
  int i=0;
  while(nh<n){
     nh=nh*2;
     i++;
  }
  return(i);
}
