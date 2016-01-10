#include <stdio.h>

double E_LJpartial(double **x, int N, int index)
{
  int i, cmpt;
  double xi[3], xij[3], rijsq;
  double sr2, sr6, Vij, V;

  V = 0.0;
  for(i=0;i<N;i++){
    if(((i-index)*(i-index))>9){
      rijsq = 0.0;
      for(cmpt=0;cmpt<3;cmpt++){
	xij[cmpt] = x[cmpt][index]-x[cmpt][i];
	rijsq += xij[cmpt]*xij[cmpt];
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
      
