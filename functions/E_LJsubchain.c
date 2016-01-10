#include <stdio.h>

double E_LJsubchain(double **x, int N, int n)
{
  double xi[3], xij[3], rijsq;
  double sr2, sr6, Vij, V;
  int i, j, k;

  V = 0.0;

  for(i=n;i<N;i++){    
      for(k=0;k<3;k++)
	xi[k] = x[k][i];
      for(j=0;j<i-3;j++){
	rijsq = 0.0;
	for(k=0;k<3;k++){
	  xij[k] = x[k][j]-xi[k];
	  rijsq += xij[k]*xij[k];
	}
	
	if(rijsq < 9.0){
	  sr2 = 1.0/rijsq;
	  sr6 = sr2*sr2*sr2;
	  Vij = sr6*(sr6-1.0);	  
	  V += Vij;
	}
      }
  }
  return 4.0*V;

}
