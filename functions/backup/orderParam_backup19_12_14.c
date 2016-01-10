/*orderParam.c
  A script that takes a set of coordinates and returns the size of the largest crystal nucleus, where crystal particles are defined using the parameter in Yi et. al. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RP2 1.2 //the larger radius for particles involved in computint p2[i]
#define RTH 1.2 //the smaller radius for the separation distance for particles in the crystal phase
#define P2THLOCAL 0.9 //threshold over which the local p2[i] must be for particle i to be considered "locally-straight"
#define P2TH 0.9 //threshold over which p2[i] must be for particle i to be considered "crystal-like"

int orderParam(double **x, int N, int SUB_CRYST, int DENSE_STATES, int output, int *particleState, int *foldCentre)
{
  int i, j, k, cmpt;
  double rsq, rijsq[N][N];
  int NTneighbours[N], NTindices[N][N], neighbourIndex;
  double vec1[3], vec2[3];
  double modVec1sq, modVec2sq, cth_ijsq;
  double p2[N];
  double Local3sum, Local5sum, Local7sum, sum2, parallel, mostParallel;
  int localStraight[N], crystal[N];
  int Ncryst=0, bestNTneighbour, bestLocalStraight = 3, bestSubCrystalState=0;
  int index;
  int absolute_zero = 1;

  for(i=0;i<N;i++){
    NTneighbours[i] = 0;
    localStraight[i]=0;
    crystal[i] = 0;
    particleState[i] = 2;
  }
 
  double RP2sq = RP2*RP2;
  //double RTHsq = RTH*RTH;

  index = 0;  
  //determine the number of neighbours of each particle i and the corresponding indices k, then store in Nb[i] and neighbours[i][k] respectively.
  for(i=3;i<N-3;i++){
    k=0;
    for(j=0;j<N;j++){
      if(((i-j)*(i-j)) > 9){	
	rsq = 0.0;
	for(cmpt=0;cmpt<3;cmpt++)
	  rsq += (x[cmpt][j]-x[cmpt][i])*(x[cmpt][j]-x[cmpt][i]);
	rijsq[i][j] = rsq;
	if(rsq < RP2sq){
	  NTneighbours[i]++;
	  NTindices[i][k] = j;	  	
	  k++;
	}	
      }
    }

    if(NTneighbours[i]>0)
      absolute_zero = 0;

  /*
  for(i=0;i<N;i++){
    if(NTneighbours[i]==0)
      printf("NTneighbours[%d]: %d\n", i, NTneighbours[i]);
    else
      printf("NTneighbours[%d]: %d (", i, NTneighbours[i]);
    for(j=0;j<NTneighbours[i];j++){
      if(j==NTneighbours[i]-1)
	printf("%d)\n", NTindices[i][j]);
      else
	printf("%d, ", NTindices[i][j]);
    }
  }
  */
    Local3sum = 0.0;
    Local5sum = 0.0;
    Local7sum = 0.0;
    modVec1sq = 0.0;
    for(cmpt=0;cmpt<3;cmpt++){   
      vec1[cmpt] = x[cmpt][i+1]-x[cmpt][i-1];
      modVec1sq += vec1[cmpt]*vec1[cmpt];
    }  
 
    for(j=i-3;j<=i+3;j++){
      if(j!=i && j>0 && j<N-1){
	modVec2sq=0.0;
	for(cmpt=0;cmpt<3;cmpt++){
	  vec2[cmpt] = x[cmpt][j+1]-x[cmpt][j-1];
	  modVec2sq += vec2[cmpt]*vec2[cmpt];
	}
	cth_ijsq = (vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]);
	cth_ijsq *= cth_ijsq/(modVec1sq*modVec2sq);
	Local7sum += 0.5*(3.0*cth_ijsq-1.0);
	if((j-i)*(j-i)<9)
	  Local5sum += 0.5*(3.0*cth_ijsq-1.0);
	if((j-i)*(j-i)<4)
	  Local3sum += 0.5*(3.0*cth_ijsq-1.0);
      }
    }
    if(Local7sum/6.0 > P2THLOCAL){
      localStraight[i]=7;
      if(NTneighbours[i] > 1){
	sum2=0.0;
	//modVec1sq = 0.0;
	//for(cmpt=0;cmpt<3;cmpt++){
	//vec1[cmpt] = x[cmpt][i+1]-x[cmpt][i-1];
	//modVec1sq += vec1[cmpt]*vec1[cmpt];
	//}
	mostParallel = 0.0;
	for(j=0;j<NTneighbours[i];j++){
	  neighbourIndex = NTindices[i][j];
	  if(neighbourIndex > 0 && neighbourIndex < N-1){
	    modVec2sq=0.0;
	    for(cmpt=0;cmpt<3;cmpt++){
	      vec2[cmpt] = x[cmpt][neighbourIndex+1]-x[cmpt][neighbourIndex-1];
	      modVec2sq += vec2[cmpt]*vec2[cmpt];
	    }
	    cth_ijsq = (vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]);
	    cth_ijsq *= cth_ijsq/(modVec1sq*modVec2sq);
	    parallel = 0.5*(3.0*cth_ijsq-1.0);
	    if(parallel > mostParallel){
	      mostParallel = parallel;
	      bestNTneighbour = NTindices[i][j];
	    }
	    //sum2 += 0.5*(3.0*cth_ijsq-1.0);
	  }
	}
	//p2[i] = sum1/6.0+sum2/(1.0*NTneighbours[i]);
	//if(p2[i] > P2TH){
	if(mostParallel>P2TH){
	  crystal[i]=1;
	  particleState[i]=4;
	  Ncryst++;
	}
      }
    }
    else{
      if(Local5sum/4.0 > P2THLOCAL)
	localStraight[i]=5;
      else{
	if(Local3sum/2.0 > P2THLOCAL)
	  localStraight[i]=3;
      }
    }
  }

  if(Ncryst==1){
    for(i=0;i<N;i++){
      if(crystal[i]==1){
	//foldCentreDist[(int)(0.5*(bestNTneighbour+i))]++;
	*foldCentre = (int)(0.5*(bestNTneighbour+i));
	break;
      }
    }
  }

  if(Ncryst>0){
    index = Ncryst+3*SUB_CRYST;
  }
  else{
    for(i=3;i<N-3;i++){
      if(localStraight[i] > bestLocalStraight)
	bestLocalStraight = localStraight[i];
    }
    for(i=3;i<N-3;i++){
      if(localStraight[i] == bestLocalStraight){
	if(NTneighbours[i]>bestSubCrystalState)
	  bestSubCrystalState = NTneighbours[i];
      }
    }
    if(bestSubCrystalState < SUB_CRYST)
      index = bestSubCrystalState;
    else
      index = SUB_CRYST;
    if(bestSubCrystalState > 0){
      switch(bestLocalStraight){
      case 7: 
	index += 2*SUB_CRYST;
	break;
      case 5: 
	index += SUB_CRYST;
	break;
      default: 
	index += 0;
	break;
      }
    }

      
    if(index > 0){
      for(i=3;i<N-3;i++){
	if(NTneighbours[i]==index)
	  particleState[i]=3;
      }
    }
  }

  if(output==1){
    printf("Ncryst: %d\n", Ncryst);
    printf("bestSubCrystalState: %d\n", bestSubCrystalState);
    printf("index: %d\n", index);
    printf("absolute_zero: %d\n", absolute_zero);
  }

  if(absolute_zero == 0)
    return index;
  else
    return -1;

}
