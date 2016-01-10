#include <math.h>
#include <stdio.h>

void crankShaft(double x[][500], int N, int n, double phi)
{
  int i, j, k;
  double U[3];
  double xi, yi, zi;
  double dx,dy,dz;
  double modU, modD;
  double transX, transY,transZ;
  double cp = cos(phi);
  double omcp = 1.0-cp;
  double sp = sin(phi);
  double rotationMatrix[3][3];
  
  for(j=0;j<3;j++)
    U[j] = x[j][n+1]-x[j][n-1];
  modU = sqrt(U[0]*U[0]+U[1]*U[1]+U[2]*U[2]);
  for(j=0;j<3;j++)
    U[j] /= modU;

  rotationMatrix[0][0] = cp + U[0]*U[0]*omcp;
  rotationMatrix[0][1] = U[0]*U[1]*omcp-U[2]*sp;
  rotationMatrix[0][2] = U[0]*U[2]*omcp+U[1]*sp;
  rotationMatrix[1][0] = U[0]*U[1]*omcp+U[2]*sp;
  rotationMatrix[1][1] = cp+U[1]*U[1]*omcp;
  rotationMatrix[1][2] = U[1]*U[2]*omcp-U[0]*sp;
  rotationMatrix[2][0] = U[0]*U[2]*omcp-U[1]*sp;
  rotationMatrix[2][1] = U[1]*U[2]*omcp+U[0]*sp;
  rotationMatrix[2][2] = cp+U[2]*U[2]*omcp;

  double r_bond[3];
  for(i=0;i<3;i++)
    r_bond[i] = x[i][n]-x[i][n-1];
  double x_rotated[3] = {0.0};

  for(i=0;i<3;i++){
    for(j=0;j<3;j++)
      x_rotated[i] += rotationMatrix[i][j]*r_bond[j];  //x[j][n];    
    x_rotated[i] += x[i][n-1];
  }
  for(i=0;i<3;i++)
    x[i][n] = x_rotated[i];
}
