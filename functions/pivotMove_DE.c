#include <math.h>
#include <stdio.h>
#include "../flexiChain.h"

double pivotMove_DE(double x[][500],double x_old[][500], int N, int pivotParticle, double lambdaSq)
{

  double E_lost=0.0;
  double E_gain=0.0;
  double rijSq;

  


    int i,j;
    double tempE;

  //compute E_lost

  for(j=0; j<= pivotParticle-1 ; j++){
    rijSq = (x_old[0][pivotParticle+1] - x_old[0][j]) * (x_old[0][pivotParticle+1] - x_old[0][j])
      +     (x_old[1][pivotParticle+1] - x_old[1][j]) * (x_old[1][pivotParticle+1] - x_old[1][j])
      +     (x_old[2][pivotParticle+1] - x_old[2][j]) * (x_old[2][pivotParticle+1] - x_old[2][j]);      

      E_lost    +=   E_swSingle( rijSq  ,  lambdaSq) ;
  }

  for(i = pivotParticle +2  ;   i <= N-1    ;   i++){
    for(j=0; j<=pivotParticle ; j++){
      rijSq = (x_old[0][i] - x_old[0][j]) * (x_old[0][i] - x_old[0][j])
	+     (x_old[1][i] - x_old[1][j]) * (x_old[1][i] - x_old[1][j])
	+     (x_old[2][i] - x_old[2][j]) * (x_old[2][i] - x_old[2][j]);
      E_lost    +=   E_swSingle( rijSq  ,  lambdaSq) ;
    }
  }
    


  

  //compute E_gain

  for(j=0; j<= pivotParticle-1 ; j++){
    rijSq = (x[0][pivotParticle+1] - x[0][j]) * (x[0][pivotParticle+1] - x[0][j])
      +     (x[1][pivotParticle+1] - x[1][j]) * (x[1][pivotParticle+1] - x[1][j])
      +     (x[2][pivotParticle+1] - x[2][j]) * (x[2][pivotParticle+1] - x[2][j]);
    /*
    tempE = E_swSingle(rijSq,lambdaSq);
    if(tempE < 0.0)
      E_gain += tempE;
    else
      return 1000.0;
      */
      E_gain    +=   E_swSingle( rijSq  ,  lambdaSq) ;
  }

  for(i = pivotParticle +2  ;   i <= N-1    ;   i++){
    for(j=0; j<=pivotParticle ; j++){
      rijSq = (x[0][i] - x[0][j]) * (x[0][i] - x[0][j])
	+     (x[1][i] - x[1][j]) * (x[1][i] - x[1][j])
	+     (x[2][i] - x[2][j]) * (x[2][i] - x[2][j]);
/*
  tempE = E_swSingle(rijSq,lambdaSq);
    if(tempE < 0.0)
      E_gain += tempE;
    else
      return 1000.0;
      */
      E_gain    +=   E_swSingle( rijSq  ,  lambdaSq) ;
    }
  }

  return E_gain - E_lost;

}

