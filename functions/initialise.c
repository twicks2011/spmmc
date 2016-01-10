#include <math.h>
#include <stdio.h>
#include <string.h>

void initialise(double x[][500], int N, char *coordsFile)
{
  int i,cmpt;  
  if(strcmp(coordsFile,"zig-zag")==0){    
     //double r = 1.1225; 
    //double r=2.0;
    double r=1.0;
    //double r = 0.382;

    x[0][0] = 0.0;
    x[1][0] = 0.0;
    x[2][0] = 0.0;
    x[0][1] = r*sin(0.955);
    x[1][1] = 0.0;
    x[2][1] = r*cos(0.955);
    for(i=2;i<N;i++){
      x[0][i] = x[0][1]*i;
      x[1][i] = 0.0;
      x[2][i] = x[2][1]*(i%2);
    }
  }
  else{        
    FILE *coordsPtr;
    coordsPtr = fopen(coordsFile,"r");    
    for(i=0;i<N;i++){
      for(cmpt=0;cmpt<3;cmpt++)
        fscanf(coordsPtr,"%lf",&x[cmpt][i]);
    }    
    return;

  }

}
