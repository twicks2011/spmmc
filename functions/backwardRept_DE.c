#include <math.h>
#include <stdio.h>
#include "../flexiChain.h"

double backwardRep_DE(double x[][500],double x_old[][500], int N, double lambdaSq)
{
  //backward reptation move, so the x_old[N] has been removed and x[0] has been created
  double E_lost=0.0;
  double E_gain=0.0;
  double rijSq;


  int i;
  double tempE;

  //compute E_lost
  for(i = 0  ;   i <= N-3   ;   i++){

    rijSq = (x_old[0][N-1] - x_old[0][i]) * (x_old[0][N-1] - x_old[0][i])
      +     (x_old[1][N-1] - x_old[1][i]) * (x_old[1][N-1] - x_old[1][i])
      +     (x_old[2][N-1] - x_old[2][i]) * (x_old[2][N-1] - x_old[2][i]);
    
     E_lost    +=   E_swSingle( rijSq  ,  lambdaSq) ;

  }
    

  //compute E_gain
  for(i = 2  ;   i <= N-1   ;   i++){

    rijSq = (x[0][0] - x[0][i]) * (x[0][0] - x[0][i])
      +     (x[1][0] - x[1][i]) * (x[1][0] - x[1][i])
      +     (x[2][0] - x[2][i]) * (x[2][0] - x[2][i]);
    /*
    tempE = E_swSingle(rijSq,lambdaSq);
    if(tempE < 0.0)
      E_gain += tempE;
    else
      return 1000.0;
*/
    E_gain    +=   E_swSingle( rijSq  ,  lambdaSq) ;

  }
  

  return E_gain - E_lost;

}

