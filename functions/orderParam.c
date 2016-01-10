/*orderParam.c
  A script that takes a set of coordinates and returns the size of the largest crystal nucleus, where crystal particles are defined using the parameter in Yi et. al. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RP2 1.2 //the larger radius for particles involved in computing p2[i]
//#define RTH 1.2 //the smaller radius for the separation distance for particles in the crystal phase
#define P2THLOCAL 0.4 //threshold over which the local p2[i] must be for particle i to be considered "locally-straight"
#define P2TH 0.4 //threshold over which p2[i] must be for particle i to be considered "crystal-like"
#define SIZE 200

int orderParam(double x[][500], int N, int output, int *particleState,double E_sw)
{  

  int i, j, k, cmpt;
  double rsq, rijsq[SIZE][SIZE], dirVec[SIZE][SIZE];
  int NTneighbours[SIZE], NTindices[SIZE][SIZE],neighbourIndex;
  double modDirVecSq[SIZE], cth_ijsq;
  double LocalSum, parallel;
  int localStraight[SIZE], crystalParticles[SIZE], crystalIndices[SIZE];
  int Ncryst=0;
  int index=0;
  int absolute_zero = 1;
  int crystCount = 0;
  int alignedCount[SIZE];
  int alignedMax = 0;
  int alignedNeighbours[SIZE][SIZE]; 
  int alignedParts = 0;

  int locallyStraightParts = 0;

/*
//END-TO-END VECTOR PARAMETER
  return floor( sqrt( (x[0][N-1]-x[0][0])*(x[0][N-1]-x[0][0])+
		      (x[1][N-1]-x[1][0])*(x[1][N-1]-x[1][0]) +
		      (x[2][N-1]-x[2][0])*(x[2][N-1]-x[2][0]) ));
          */

  for(i=0;i<N;i++){
    NTneighbours[i] = 0;
    localStraight[i]=0;
    crystalParticles[i] = 0;
    particleState[i] = 2; //for VMD
    alignedCount[i] = 0;
    for(j=0;j<N;j++)
      alignedNeighbours[i][j] = 0;
  }
  
  //compute the direction vector for each particle (this may be optimised later)
  for(i=1;i<N-1;i++){
    modDirVecSq[i] = 0.0;
    for(cmpt=0;cmpt<3;cmpt++){
      dirVec[cmpt][i] = x[cmpt][i+1]-x[cmpt][i-1];
      modDirVecSq[i] += dirVec[cmpt][i]*dirVec[cmpt][i];
    }
  }

  double RP2sq = RP2*RP2;
  for(i=2;i<N-2;i++){
    //determine th number of non-trivial neighbours of each particle i and the corresponding indices
    k=0;
    for(j=0;j<N;j++){
      if(((i-j)*(i-j))>9){        
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
    /*
    printf("NTindices[%d]: %d", i, NTindices[i][0]);
    for(j=1;j<NTneighbours[i];j++)
      printf(", %d",NTindices[i][j]);
    printf("\n");
    */
    if(NTneighbours[i]>0){
      //Determine whether each partical is "locally straight" or not.
      LocalSum = 0.0;
      for(j=i-1;j<=i+1;j++){
      //for(j=i-2;j<=i+2;j++){			
      //for(j=i-3;j<=i+3;j++){      
        if(i!=j){          
          cth_ijsq = dirVec[0][i]*dirVec[0][j]+dirVec[1][i]*dirVec[1][j]+dirVec[2][i]*dirVec[2][j];
	        cth_ijsq *= cth_ijsq/(modDirVecSq[i]*modDirVecSq[j]);
	        LocalSum += 0.5*(3.0*cth_ijsq-1.0);
        }
      }
      //printf("%f\n",LocalSum);
      if(LocalSum/2.0>P2THLOCAL){
      //if(LocalSum/4.0>P2THLOCAL){
      //if(LocalSum/6.0>P2THLOCAL){
	      localStraight[i]=1;
	      particleState[i]=3;
	      locallyStraightParts++;
      }
    }
  }

  int parallelCount;
  for(i=2;i<N-2;i++){
    if(localStraight[i] == 1){
      parallel = 0.0;
      parallelCount = 0;
      for(j=0;j<NTneighbours[i];j++){
        neighbourIndex = NTindices[i][j];
        if(localStraight[neighbourIndex] == 1){ //***This may be too strict
          parallelCount++;
          cth_ijsq = dirVec[0][i]*dirVec[0][neighbourIndex]+dirVec[1][i]*dirVec[1][neighbourIndex]+dirVec[2][i]*dirVec[2][neighbourIndex];
          cth_ijsq *= cth_ijsq/(modDirVecSq[i]*modDirVecSq[neighbourIndex]);
          parallel += 0.5*(3.0*cth_ijsq-1.0);
        }        
      }
      if(parallelCount>0 && (parallel/(1.0*parallelCount)) > P2TH){
        crystalParticles[i] = 1;
        particleState[i] = 1; //for VMD
      }
    }
  }

  j=0;
  for(i=0;i<N;i++){
    if(crystalParticles[i] == 1){
      crystCount++;
      crystalIndices[j] = i;
      j++;
    }
  }

  for(i=0;i<crystCount;i++){
    for(j=0;j<crystCount;j++){
      if(i!=j){
        cth_ijsq = dirVec[0][crystalIndices[i]]*dirVec[0][crystalIndices[j]]+dirVec[1][crystalIndices[i]]*dirVec[1][crystalIndices[j]]+dirVec[2][crystalIndices[i]]*dirVec[2][crystalIndices[j]];
        cth_ijsq *= cth_ijsq/(modDirVecSq[crystalIndices[i]]*modDirVecSq[crystalIndices[j]]);
        parallel = 0.5*(3.0*cth_ijsq-1.0);
        if(parallel>P2TH){
          alignedNeighbours[i][crystalIndices[j]] = 1;
          alignedCount[i]++;
        }
      }
    }
    if(alignedCount[i] > alignedMax)
      alignedMax = alignedCount[i];
  }

  for(i=0;i<crystCount;i++){
    alignedParts = 0;
    if(alignedCount[i] == alignedMax){
      particleState[crystalIndices[i]] = 4;
      alignedParts++;
      for(j=0;j<N;j++){
        if(alignedNeighbours[i][j] == 1){
          particleState[j] = 4;
          alignedParts++;
        }
      }
      if(alignedParts == 1)
        particleState[crystalIndices[i]] = 1;
    }
  }

  //return locallyStraightParts;
  //return crystCount;

  if(alignedMax == 0)
    return 0;
  else
    return alignedMax +1;


}
