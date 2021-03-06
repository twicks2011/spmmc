/*orderParam.c
  A script that takes a set of coordinates and returns the size of the largest crystal nucleus, where crystal particles are defined using the parameter in Yi et. al. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RP2 1.2 //the larger radius for particles involved in computing p2[i]
#define RTH 1.2 //the smaller radius for the separation distance for particles in the crystal phase
#define P2THLOCAL 0.98 //threshold over which the local p2[i] must be for particle i to be considered "locally-straight"
#define P2TH 0.98 //threshold over which p2[i] must be for particle i to be considered "crystal-like"
#define SIZE 200

int orderParam(double x[][500], int N, int output, int *particleState)
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
  //double RTHsq = RTH*RTH;

  //index = 0;  
  //determine the number of neighbours of each particle i and the corresponding indices k, then store in Nb[i] and neighbours[i][k] respectively.
  for(i=2;i<N-2;i++){
    k=0;
    for(j=0;j<N;j++){
      if(((i-j)*(i-j)) > 9){	
	//if(((i-j)*(i-j)) > 1){	
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
    
    if(NTneighbours[i]>0){
      absolute_zero = 0;
      /*
	if(output==1){
	for(i=0;i<N;i++){
	for(j=0;j<N;j++)
	printf("%.2f ",rijsq[i][j]);
	printf("\n");
	}
	}
      */	      

      //Determine whether each partical is "locally straight" or not.
      LocalSum = 0.0;
      //for(j=i;j<=i+1;j++){			
	//for(j=i-1;j<=i+1;j++){			
      for(j=i-2;j<=i+2;j++){
      //for(j=i-3;j<=i+3;j++){
	cth_ijsq = dirVec[0][i]*dirVec[0][j]+dirVec[1][i]*dirVec[1][j]+dirVec[2][i]*dirVec[2][j];
	cth_ijsq *= cth_ijsq/(modDirVecSq[i]*modDirVecSq[j]);
	LocalSum += 0.5*(3.0*cth_ijsq-1.0);
      }
      //if(LocalSum/2.0 > P2THLOCAL){
      //if(LocalSum/1.0 > P2THLOCAL){
      if(LocalSum/4.0>P2THLOCAL){
      //if(LocalSum/6.0>P2THLOCAL){
	localStraight[i]=1;
	particleState[i]=3;
	locallyStraightParts++;
      }
    }

    //printf("localStraight[%d]: %d\n", i, localStraight[i]);

  }
  int parallelCount;
  for(i=2;i<N-2;i++){
    if(localStraight[i] == 1){
      parallel = 0.0;
      parallelCount=0;
      for(j=0;j<NTneighbours[i];j++){
	neighbourIndex = NTindices[i][j];
	if(localStraight[neighbourIndex]==1){
    parallelCount++;
	  cth_ijsq = dirVec[0][i]*dirVec[0][neighbourIndex]+dirVec[1][i]*dirVec[1][neighbourIndex]+dirVec[2][i]*dirVec[2][neighbourIndex];
	  cth_ijsq *= cth_ijsq/(modDirVecSq[i]*modDirVecSq[neighbourIndex]);
	  parallel += 0.5*(3.0*cth_ijsq-1.0);
	}
      }
      if(parallelCount>0 && (parallel/(1.0*parallelCount))>P2TH){
	crystalParticles[i] = 1;
	particleState[i] = 1; //for VMD		
      }
    }
    //printf("crystalParticles[%d]: %d\n", i, crystalParticles[i]);
  }
	  
  j=0;
  for(i=0;i<N;i++){
    if(crystalParticles[i]==1){
      crystCount++;
      crystalIndices[j]=i;
      j++;
    }
  }

  /*
  printf("crystCount: %d (",crystCount);
  for(i=0;i<crystCount;i++)
    printf("%d, ",crystalIndices[i]);
  printf(")\n");
  */

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

  //if(output == 1)
  //printf("alignedMax: %d\n", alignedMax);
  for(i=0;i<crystCount;i++){
    alignedParts = 0;
    //if(output == 1)
    //printf("alignedCount[%d]: %d\n", crystalIndices[i], alignedCount[i]);
    if(alignedCount[i]==alignedMax){
      //index++;
      particleState[crystalIndices[i]] = 4;
      alignedParts++;
      for(j=0;j<N;j++){
	if(alignedNeighbours[i][j] == 1){
	  //if(output == 1)
	  //printf("alignedNeighbours[%d][%d]: %d\n", crystalIndices[i], j, alignedNeighbours[i][j]);	
	  //index++;
	  particleState[j]=4; //for VMD
	  alignedParts++;
	}	
      }
      if(alignedParts == 1)
	particleState[crystalIndices[i]] = 1;      
      //if(output==1)
      //printf("Should stop now\n");
      break;
    }
  }
  
  if(output==1){  
    printf("crystCount: %d\n", crystCount);
    printf("alignedMax: %d\n", alignedMax);
    printf("absolute_zero: %d\n", absolute_zero);
    //for(i=0;i<N;i++)
    //printf("particleState[%d]: %d\n", i, particleState[i]);
  }
  
  //if(absolute_zero == 0){

  
  if(alignedMax==0)
    return 0;
  else
    return alignedMax+1;
  
  
  //return locallyStraightParts;

    //return index;
    //}
  //else
  //return -1;

}
