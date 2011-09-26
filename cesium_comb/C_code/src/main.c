#include <stdio.h>
#include <matrixMul.h>

int main(int argc, char *argv[])
{
  printf ("Cublas Testing programe ... \n \n");

  deviceVerify();
  runTest();          
  
  //complexTest();
  
  return 0;
  
}
