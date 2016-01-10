#include <stdio.h>

double E_LJnonCrystAv(double x[][500], int N, int *particleState)
{
  int i, j, k;
  double xi[3], xij[3], rijsq;
  double sr2, sr6, Vij, V;
  int Ncryst = 0;
  //static int computed[SIZE][SIZE];
  //int computed[N][N];
  /*
  for(i=0;i<N;i++){
    for(j=0;j<N;j++)
      computed[i][j] = 0;
  }
  */
  V = 0.0;

  for(i=0;i<N;i++){
    if(particleState[i] != 4){
      Ncryst++;
      for(k=0;k<3;k++)
	xi[k] = x[k][i];
      for(j=0;j<N;j++){
	//if(computed[i][j] == 0){ // && particleState[j]!=4){
	  //if(((i-j)*(i-j))>1){
	  //if(((i-j)*(i-j))>9){
	if(i != j ){
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
	      //computed[i][j] = 1;
	      //computed[j][i] = 1;
	}
	      //}
	  //}
      }
    }
  }
  //printf("Ncryst: %d\n", Ncryst);
  return 4.0*V;
}
