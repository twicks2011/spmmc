#include <stdio.h>
#include <math.h>

#define PI 3.141592653589793238462643383279502884197

void positions2bonds(double x[][500], double *theta, double *phi, int N)
{
  double v1x,v1y,v1z,v2x,v2y,v2z;
  double modV2, thetaCompl,ctheta, cphi;
  int i, j;
  double r;

  for(i=1;i<N;i++){
    v1x = x[0][i]-x[0][i-1];
    v1y = x[1][i]-x[1][i-1];
    v1z = x[2][i]-x[2][i-1];

    r = 1.0;//qrt(v1x*v1x + v1y*v1y + v1z*v1z);       
    //if(i==0)
    //theta[i] = acos(v1z/r[i]);
    //else{
    if(i<N-1){
      v2x = x[0][i+1]-x[0][i];
      v2y = x[1][i+1]-x[1][i];
      v2z = x[2][i+1]-x[2][i];
      modV2 = 1.0; //sqrt(v2x*v2x+v2y*v2y+v2z*v2z);
       
      ctheta = (v1x*v2x+v1y*v2y+v1z*v2z)/(r*modV2);       
      theta[i] = PI -  acos(ctheta);
    }
  }

  double rij[3], rjk[3], rkl[3], nijk[3], njkl[3], modN1, modN2;  

  for(i=0;i<N-1;i++){
    if(i<2)
      phi[i] = 0.0;
    else{
      for(j=0;j<3;j++){
	rij[j] = x[j][i-1]-x[j][i-2];
	rjk[j] = x[j][i]-x[j][i-1];
	rkl[j] = x[j][i+1]-x[j][i];
      }

      //printf("rij[%d] = (%.3f,%.3f,%.3f)\n",i,rij[0],rij[1],rij[2]);
      //printf("rjk[%d] = (%.3f,%.3f,%.3f)\n",i,rjk[0],rjk[1],rjk[2]);
      //printf("rkl[%d] = (%.3f,%.3f,%.3f)\n",i,rkl[0],rkl[1],rkl[2]);

      nijk[0] = rij[1]*rjk[2]-rij[2]*rjk[1];
      nijk[1] = rij[2]*rjk[0]-rij[0]*rjk[2];
      nijk[2] = rij[0]*rjk[1]-rij[1]*rjk[0];

      njkl[0] = rjk[2]*rkl[1]-rjk[1]*rkl[2];
      njkl[1] = rjk[0]*rkl[2]-rjk[2]*rkl[0];
      njkl[2] = rjk[1]*rkl[0]-rjk[0]*rkl[1];             

      //  printf("nijk[%d] = (%.3f,%.3f,%.3f)\n",i,nijk[0],nijk[1],nijk[2]);
      //  printf("njkl[%d] = (%.3f,%.3f,%.3f)\n",i,njkl[0],njkl[1],njkl[2]);

      modN1 = sqrt(nijk[0]*nijk[0]+nijk[1]*nijk[1]+nijk[2]*nijk[2]);
      modN2 = sqrt(njkl[0]*njkl[0]+njkl[1]*njkl[1]+njkl[2]*njkl[2]);

      cphi = (nijk[0]*njkl[0]+nijk[1]*njkl[1]+nijk[2]*njkl[2])/(modN1*modN2);
      //printf("cphi: %f\n", cphi);
      if(cphi > 1.0)
	cphi = 1.0;
      else{
	if(cphi < -1.0)
	  cphi = -1.0;
      }
      phi[i] = acos(cphi);
	
    }
  }
}
