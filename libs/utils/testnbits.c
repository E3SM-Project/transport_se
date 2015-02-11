#include <stdio.h>
int main(){

int nbits(int);

int nb;

int n;
for(n=0;n<=100;n++){
   nb=nbits(n);
   printf("n=%d nbits=%d\n",n,nb);
   }

}
