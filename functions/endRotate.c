#include <math.h>
#include <stdio.h>

void endRotate(double x[][500],int N, int end, double rndnum_1, double rndnum_2, double theta)
{
	int i, cmpt;

	double rij[3], rjk[3], nijk[3];
  double rndnum_1sq = rndnum_1*rndnum_1;
  double rndnum_2sq = rndnum_2*rndnum_2;
/*
  if(end==0){
  	 for(cmpt=0;cmpt<3;cmpt++){
	   	 rij[cmpt] = x[cmpt][1]-x[cmpt][0];
	 	   rjk[cmpt] = x[cmpt][2]-x[cmpt][1];
	   }
  }
  else{
    for(cmpt=0;cmpt<3;cmpt++){
      rij[cmpt] = x[cmpt][N-2]-x[cmpt][N-3];
      rjk[cmpt] = x[cmpt][N-1]-x[cmpt][N-2];
    }
  }

	nijk[0] = rij[1]*rjk[2]-rij[2]*rjk[1];
  nijk[1] = rij[2]*rjk[0]-rij[0]*rjk[2];
  nijk[2] = rij[0]*rjk[1]-rij[1]*rjk[0];	

    double modNijk = sqrt(nijk[0]*nijk[0]+nijk[1]*nijk[1]+nijk[2]*nijk[2]);
    for(cmpt=0;cmpt<3;cmpt++)
    	nijk[cmpt] /= modNijk;
*/

    if(end == 0){
      for(cmpt=0;cmpt<3;cmpt++)
        rij[cmpt] = x[cmpt][0]-x[cmpt][1];
    }
    else{
      for(cmpt=0;cmpt<3;cmpt++)
        rij[cmpt] = x[cmpt][N-1]-x[cmpt][N-2];
    }

    double randPt[3];
    randPt[0] = 2.0*rndnum_1*sqrt(1.0-rndnum_1sq-rndnum_2sq);
    randPt[1] = 2.0*rndnum_2*sqrt(1.0-rndnum_1sq-rndnum_2sq);
    randPt[2] = 1.0-2.0*(rndnum_1sq+rndnum_2sq);

    for(cmpt=0;cmpt<3;cmpt++){
      if(end==0)
        x[cmpt][0] = x[cmpt][1]+randPt[cmpt];
      else
        x[cmpt][N-1] = x[cmpt][N-2]+randPt[cmpt];
    }


    double rotationMatrix[3][3];
    double cp = cos(theta);
    double omcp = 1.0-cp;
    double sp = sin(theta);

    rotationMatrix[0][0] = cp + randPt[0]*randPt[0]*omcp;
  	rotationMatrix[0][1] = randPt[0]*randPt[1]*omcp-randPt[2]*sp;
  	rotationMatrix[0][2] = randPt[0]*randPt[2]*omcp+randPt[1]*sp;
  	rotationMatrix[1][0] = randPt[0]*randPt[1]*omcp+randPt[2]*sp;
  	rotationMatrix[1][1] = cp+randPt[1]*randPt[1]*omcp;
	  rotationMatrix[1][2] = randPt[1]*randPt[2]*omcp-randPt[0]*sp;
  	rotationMatrix[2][0] = randPt[0]*randPt[2]*omcp-randPt[1]*sp;
  	rotationMatrix[2][1] = randPt[1]*randPt[2]*omcp+randPt[0]*sp;
  	rotationMatrix[2][2] = cp+randPt[2]*randPt[2]*omcp;

    double rotated_bond[3]={0.0};
    double modR=0.0;

    if(end == 0){    	
    	//modR = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
    	for(cmpt=0;cmpt<3;cmpt++){    		
    		for(i=0;i<3;i++)
    			rotated_bond[cmpt] += rotationMatrix[cmpt][i]*rij[i];
    		x[cmpt][0] = x[cmpt][1] + rotated_bond[cmpt];
    	}

    }
    else{    	
      //modR = sqrt(rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2]);
    	for(cmpt=0;cmpt<3;cmpt++){
    		for(i=0;i<3;i++)
    			rotated_bond[cmpt] += rotationMatrix[cmpt][i]*rij[i];		
    		x[cmpt][N-1] = x[cmpt][N-2] + rotated_bond[cmpt];
    	}
    }       
}