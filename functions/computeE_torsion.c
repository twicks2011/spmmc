#include <math.h>
#include <stdio.h>

double computeE_torsion(double *torsionAngles, int N, int output, int state)
{
  int i;
  double  V = 0.0;
  double k1 = 14.264, k2 = -7.74, k3 = 28.913;
  //double k1 = 12.9, k2 = -10.2, k3 = 11.1;
  //double k1 = 0.0, k2 = 0.0, k3 = 20.0;
 
  /*
  FILE *stateE_TorPtr;
  if(output == 1){
    char stateE_TorFile[100];
    int dummyInt = sprintf(stateE_TorFile,"stateE_Tor/state%d.dat", state);
    stateE_TorPtr = fopen(stateE_TorFile, "w");
  }
  */
  //for(i=0;i<N;i++)
  //printf("phiINSIDE[%d]: %f\n",i,torsionAngles[i]); 

  for(i=2;i<N;i++){
    //if(output == 1)
    //fprintf(stateE_TorPtr, "%d %f\n", i, 0.5*(k1*(1.0-cos(torsionAngles[i])) + k2*(1.0-cos(2*torsionAngles[i])) + k3*(1.0-cos(3*torsionAngles[i]))));
    V += k1*(1.0-cos(torsionAngles[i])) + k2*(1.0-cos(2*torsionAngles[i])) + k3*(1.0-cos(3*torsionAngles[i]));          
  }
  //if(output == 1)
  //fclose(stateE_TorPtr);

  return 0.5*V;
}
