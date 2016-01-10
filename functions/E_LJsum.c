#include <stdio.h>

double E_LJsum(double x[][500], int N)
{
  int i,j,k;
  double xi[3], xij[3], rijsq;
  double sr2, sr6, Vij, V;

  /*
  int computed[N][N];
  for(i=0;i<N;i++){
    for(j=0;j<N;j++)
      computed[i][j] = 0;
  }
  */

  V = 0.0;

  //for(i=0;i<N-4;i++){
  for(i=0;i<N;i++){
    for(k=0;k<3;k++)
      xi[k] = x[k][i];
    //for(j=i+4;j<N;j++){
    for(j=i+1;j<N;j++){
      rijsq = 0.0;
      for(k=0;k<3;k++){
	xij[k] = x[k][j]-xi[k];
	rijsq += xij[k]*xij[k];
      }
      
      //if(rijsq < 9.0){     
	sr2 = 1.0/rijsq;
	sr6 = sr2*sr2*sr2;
	Vij = sr6*(sr6-1.0);
	V += Vij;
	//}
    }
  }
  return 4.0*V;
}
