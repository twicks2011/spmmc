/*orderParam.c
  A script that takes a set of coordinates and returns the size of the largest crystal nucleus, where crystal particles are defined using the parameter in Yi et. al. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RP2 1.2 //the larger radius for particles involved in computing p2[i]
#define RTH 1.2 //the smaller radius for the separation distance for particles in the crystal phase
#define P2THLOCAL 0.7 //threshold over which the local p2[i] must be for particle i to be considered "locally-straight"
#define P2TH 0.7 //threshold over which p2[i] must be for particle i to be considered "crystal-like"

int orderParam(double **x, int N, int SUB_CRYST, int DENSE_STATES, int output, int *particleState, int *foldCentre)
{
  int i, j, k, cmpt;
  double rsq, rijsq[N][N], dirVec[N][N];
  int NTneighbours[N], NTindices[N][N], parallelNeighbours[N][N],neighbourIndex;
  //double vec1[3], vec2[3];
  //double modVec1sq, modVec2sq, cth_ijsq;
  double modDirVecSq[N], cth_ijsq;
  double p2[N];
  double Local3sum, Local5sum, Local7sum, sum2, parallel, mostParallel;
  int localStraight[N], crystal[N];
  int Ncryst=0, bestNTneighbour[N], bestParticle=0, firstComponent[N], secondComponent[N], dummyCmpt[N], cmptSum, cmptDiff, cmptSumMax, cmptDiffMin, bestSubCrystalState=0;
  int index;
  int absolute_zero = 1;

  int crystCount = 0;

  for(i=0;i<N;i++){
    NTneighbours[i] = 0;
    localStraight[i]=0;
    firstComponent[i]=0;
    secondComponent[i]=0;
    dummyCmpt[i]=0;
    crystal[i] = 0;
    particleState[i] = 2;
  }

  for(i=1;i<N-1;i++){
    modDirVecSq[i] = 0.0;
    for(cmpt=0;cmpt<3;cmpt++){
      dirVec[cmpt][i] = x[cmpt][i+1]-x[cmpt][i-1];
      modDirVecSq[i] += dirVec[cmpt][i]*dirVec[cmpt][i];
    }
  }

  double RP2sq = RP2*RP2;
  //double RTHsq = RTH*RTH;

  index = 0;  
  //determine the number of neighbours of each particle i and the corresponding indices k, then store in Nb[i] and neighbours[i][k] respectively.
  for(i=1;i<N-1;i++){
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
    if(output==1){
      for(i=0;i<N;i++){
	for(j=0;j<N;j++)
	  printf("%.2f ",rijsq[i][j]);
	printf("\n");
      }
    }
    */	      

    //Determine the "local straightness" of each particle. Each particle i will be given an index [0,3,5,7], indicating the length of the straight segment, with particle i at the centre.
    Local3sum = 0.0;
    Local5sum = 0.0;
    Local7sum = 0.0;
    //modVec1sq = 0.0; 
    for(j=i-3;j<=i+3;j++){
      if(j!=i && j>0 && j<N-1){
	//modVec2sq=0.0;			
	cth_ijsq = dirVec[0][i]*dirVec[0][j]+dirVec[1][i]*dirVec[1][j]+dirVec[2][i]*dirVec[2][j];
	cth_ijsq *= cth_ijsq/(modDirVecSq[i]*modDirVecSq[j]);
	Local7sum += 0.5*(3.0*cth_ijsq-1.0);
	if((j-i)*(j-i)<9)
	  Local5sum += 0.5*(3.0*cth_ijsq-1.0);
	if((j-i)*(j-i)<4)
	  Local3sum += 0.5*(3.0*cth_ijsq-1.0);
      }
    }
    if(Local7sum/6.0 > P2THLOCAL){
      localStraight[i]=7;
      //if(bestLocalStraight<7)
      //bestLocalStraight=7;
    }
    else{
      if(Local5sum/4.0 > P2THLOCAL){
	localStraight[i]=5;
	//if(bestLocalStraight<5)
	//bestLocalStraight=5;
      }
      else{
	if(Local3sum/2.0 > P2THLOCAL){
	  localStraight[i]=3;
	  //if(bestLocalStraight<3)
	  //bestLocalStraight=3;
	}
      }
    }
  }

  for(i=1;i<N-1;i++){
    //if(localStraight[i]==bestLocalStraight && NTneighbours[i]>0){
    //if(localStraight[i]>0 && NTneighbours[i]>0){
    if(NTneighbours[i]>0){
      if(localStraight[i] == 0)
	firstComponent[i] = 1;
      else
	firstComponent[i] = localStraight[i];     
      //modVec1sq = 0.0;
      //for(cmpt=0;cmpt<3;cmpt++){   
      //vec1[cmpt] = x[cmpt][i+1]-x[cmpt][i-1];
      //modVec1sq += vec1[cmpt]*vec1[cmpt];
      //}
      //mostParallel = 0.0;
      cmptSumMax = 0;
      cmptDiffMin = 50;
      for(j=0;j<NTneighbours[i];j++){
	neighbourIndex = NTindices[i][j];
	if(neighbourIndex > 0 && neighbourIndex < N-1 && ((i-neighbourIndex)*(i-neighbourIndex))>25){
	  //modVec2sq=0.0;
	  //for(cmpt=0;cmpt<3;cmpt++){
	  //vec2[cmpt] = x[cmpt][neighbourIndex+1]-x[cmpt][neighbourIndex-1];
	  //modVec2sq += vec2[cmpt]*vec2[cmpt];
	  //}
	  cth_ijsq = dirVec[0][i]*dirVec[0][neighbourIndex]+dirVec[1][i]*dirVec[1][neighbourIndex]+dirVec[2][i]*dirVec[2][neighbourIndex];
	  cth_ijsq *= cth_ijsq/(modDirVecSq[i]*modDirVecSq[neighbourIndex]);
	  //cth_ijsq = (vec1[0]*vec2[0]+vec1[1]*vec2[1]+vec1[2]*vec2[2]);
	  //cth_ijsq *= cth_ijsq/(modVec1sq*modVec2sq);
	  parallel = 0.5*(3.0*cth_ijsq-1.0);
	  //if(parallel>mostParallel)
	  //mostParallel = parallel;
	  //if(parallel>P2TH && localStraight[NTindices[i][j]]>bestLocalStraightNeighbour){
	  if(parallel>P2TH){
	    //if(output==1)
	    //printf("parallel[%d][%d]: %f\n", i, neighbourIndex, parallel);
	    parallelNeighbours[i][j]=1;
	    if(localStraight[neighbourIndex]==0)
	      dummyCmpt[j] = 1;
	    else
	      dummyCmpt[j] = localStraight[neighbourIndex];
	  }	  
	  else {
	    parallelNeighbours[i][j]=0;
	    dummyCmpt[j] = 0;
	  }
	  cmptSum = firstComponent[i]+dummyCmpt[j];
	  if(cmptSum>cmptSumMax)
	    cmptSumMax = cmptSum;
	}
      }     
      for(j=0;j<NTneighbours[i];j++){
	if(firstComponent[i]+dummyCmpt[j]==cmptSumMax){
	  cmptDiff = (firstComponent[i]-dummyCmpt[j])*(firstComponent[i]-dummyCmpt[j]);
	  if(cmptDiff<cmptDiffMin){
	    secondComponent[i] = dummyCmpt[j];
	    bestNTneighbour[i]=NTindices[i][j];
	    cmptDiffMin = cmptDiff;
	  }
	}
      }
    }       
  }
  /* 
 if(output==1){
    for(i=0;i<N;i++){
      if(NTneighbours[i]==0)
	printf("NTneighbours[%d]: %d\n", i, NTneighbours[i]);
      else
	printf("NTneighbours[%d]: %d (", i, NTneighbours[i]);
      for(j=0;j<NTneighbours[i];j++){
	if(j==NTneighbours[i]-1){
	  if(parallelNeighbours[i][j]==1)
	    printf("%dp)\n", NTindices[i][j]);
	  else
	    printf("%d)\n", NTindices[i][j]);
	}
	else{
	  if(parallelNeighbours[i][j]==1)
	    printf("%dp, ", NTindices[i][j]);
	  else	    
	    printf("%d, ", NTindices[i][j]);
	}
      }
    }
    for(i=0;i<N;i++)
      printf("localStraight[%d]: %d\n", i, localStraight[i]);
  }
  */	  

  cmptDiffMin = 50;
  cmptSumMax = 0;
  for(i=1;i<N-1;i++){
    if(NTneighbours[i]>0){
      if(secondComponent[i]>firstComponent[i]){
	dummyCmpt[i] = firstComponent[i];
	firstComponent[i] = secondComponent[i];
	secondComponent[i] = dummyCmpt[i];
      }
    }
    if(firstComponent[i]+secondComponent[i]>cmptSumMax)
      cmptSumMax = firstComponent[i]+secondComponent[i];
  }
  for(i=1;i<N-1;i++){
    if(firstComponent[i]+secondComponent[i]==cmptSumMax){
      if(firstComponent[i]-secondComponent[i]<cmptDiffMin){
	bestParticle = i;
	cmptDiffMin = firstComponent[i]-secondComponent[i];
      } 
      //if(firstComponent[i]>=3 && secondComponent[i]>=3){
      if(SUB_CRYST>0){
	if(firstComponent[i]==7 && secondComponent[i]==7){
	  particleState[i] = 4;
	  crystCount++;
	}
      }
      else{
	if(firstComponent[i]>=3 && secondComponent[i]>=3){
	  particleState[i] = 4;
	  crystCount++;
	}
      }
    }
  }

  if(SUB_CRYST>0){
    int cOne = firstComponent[bestParticle], cTwo = secondComponent[bestParticle];
    if(crystCount==0){
      if(cOne==3 && cTwo == 3)
	index = 1;
      else{
	if(cOne == 5){
	  if(cTwo == 3) 
	    index = 2;
	  else{
	    if(cTwo == 5)
	      index = 3;
	  }
	}
	else{
	  if(cOne == 7){
	    if(cTwo == 3)
	      index = 4;
	    else{
	      if(cTwo == 5)
		index = 5;
	    }
	  }
	  else
	    index = 0;
	}
      }
    }
    else
      index = 5 + crystCount;
  }
  else
    index = crystCount;
	  
  



  /*
  if(output==1){
    for(i=0;i<N;i++)
      printf("Pair-state[%d]: (%d,%d)\n",i,firstComponent[i],secondComponent[i]);
    printf("bestParticle: %d\n", bestParticle);
    printf("bestNTneighbour: %d\n", bestNTneighbour[bestParticle]);
  }
  */
 

  //particleState[bestParticle] = 4;
  //particleState[bestNTneighbour[bestParticle]]=3;
  //index = 10*firstComponent[bestParticle]+secondComponent[bestParticle];
  //index = crystCount;

  /*
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
  */
  /*
  if(Ncryst==1){
    for(i=0;i<N;i++){
      if(crystal[i]==1){
	//foldCentreDist[(int)(0.5*(bestNTneighbour+i))]++;
	*foldCentre = (int)(0.5*(bestNTneighbour+i));
	break;
      }
    }
  }
  */
  /*
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
  */
  
  if(output==1){
    //printf("Ncryst: %d\n", Ncryst);
    printf("crystCount: %d\n", crystCount);
    //printf("bestSubCrystalState: %d\n", bestSubCrystalState);
    printf("index: %d\n", index);
    printf("absolute_zero: %d\n", absolute_zero);
  }
  
  if(absolute_zero == 0)
    return index;
  else
    return -1;

}
