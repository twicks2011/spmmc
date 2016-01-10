#include <math.h>
#include <stdio.h>

double E_swSingle(double rijsq,  double lambdaSq)
{

  
  if (rijsq < 1) return 1000.0;
  if (rijsq > lambdaSq) return 0;
  return -1.0;

}
