/*main.c 
Single chain model with Parallel Tempering implemented using OpenMP. A range of temperatures is defined at the start of the simulation, then we use multithreading to run parallel simulations at each temperature, swapping information between chains at regular intervals.
*/

/* PREPROCESSOR */

//headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "flexiChain.h"
#include <omp.h>


//constants
#define K_B 1.3806488e-23
#define AVOGADRO_CNST 6.0221415e23
#define EPSILON 469
#define PI 3.141592653589793238462643383279502884197
#define PHI_BIG 2.094395102
#define THRESHOLD (int)1e4
#define N_Cstar 4

int main(int argc, char *argv[])
{
  int N;
  double LAMBDA, LAMBDA_SQ;
  int TEMP_INTEREST;
  int NUM_TEMPS;
  //double RAT;
  char tempsFile[100];
  int fileInput;
  char configsInput[100];
  char moveSizesDir[100];
  double PHI_MAX;
  double ALPHA;
  double BETA;
  double GAMMA;
  double sigma;
  int stretchStart;
  int stretchFinish;
  int Nmin;
  int Nmax;
  double kappa;
  double kappaRHS;
  //char updateInfoDir[100];
  int SUB_CRYST;
  int update;
  double tolerance;
  double PRINT_INTERVAL;
  double PHI_UPDATE_INTERVAL;
  double UPDATE_INTERVAL;
  double SWAP_INTERVAL;
  int NUM_SWAPS;
  double NUM_SWAP_INTERVALS;
  double RESET_INTERVAL;
  double SMALL_ANGLE_MOVES;
  double REPTATION_MOVES;
  double CRANK_MOVES;
  double END_ROT_MOVES;
  double END_BRIDGE_MOVES;
  char outputDir[100];
  char weightFile[100];
  char VMDdir[100];
  char initCoords[100];  
  

  FILE *inputPtr;
  char buffer[100];

  if(argc<2){
    printf("Usage: ./singleChainOMP [input filename]\n");
    exit(EXIT_FAILURE);
  }
  
  if((inputPtr = fopen(argv[1],"r")) == NULL){
    printf("Cannot open file %s\n",argv[1]);
    exit(EXIT_FAILURE);
  }

  printf("Reading from the input file: %s\n", argv[1]);
  fscanf(inputPtr, "%d\n", &N);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &LAMBDA);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%d\n", &TEMP_INTEREST);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%d\n", &NUM_TEMPS);
  fgets(buffer, 100, inputPtr);
  //fscanf(inputPtr, "%lf\n", &RAT);
  fscanf(inputPtr, "%s\n", tempsFile);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%d\n", &fileInput);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%s\n", configsInput);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%s\n",moveSizesDir);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &PHI_MAX);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &ALPHA);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &BETA);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &GAMMA);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &sigma);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%d\n", &stretchStart);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%d\n", &stretchFinish);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%d\n", &Nmin);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%d\n", &Nmax);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &kappa);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &kappaRHS);
  fgets(buffer, 100, inputPtr);
  //fscanf(inputPtr, "%s\n", updateInfoDir);
  //fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%d\n", &update);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &tolerance);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%le\n", &PRINT_INTERVAL);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%le\n", &PHI_UPDATE_INTERVAL);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%le\n", &UPDATE_INTERVAL);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%le\n", &SWAP_INTERVAL);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%d\n", &NUM_SWAPS);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%le\n", &NUM_SWAP_INTERVALS);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%le\n", &RESET_INTERVAL);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &SMALL_ANGLE_MOVES);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &REPTATION_MOVES);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &CRANK_MOVES);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &END_ROT_MOVES);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%lf\n", &END_BRIDGE_MOVES);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%s\n",outputDir);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%s\n",weightFile);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%s\n",VMDdir);
  fgets(buffer, 100, inputPtr);
  fscanf(inputPtr, "%s\n",initCoords);
  fgets(buffer, 100, inputPtr);
  fclose(inputPtr);
/*
  if(Nmax > N){
    printf("ERROR! Nmax cannot exceed N\n");
    exit(EXIT_FAILURE);
  }
*/

  LAMBDA_SQ = LAMBDA * LAMBDA;

  
  double propSum = SMALL_ANGLE_MOVES+REPTATION_MOVES+CRANK_MOVES+END_ROT_MOVES+END_BRIDGE_MOVES;
  int propError = (propSum-1.0)*(propSum-1.0) < 1e-8 ? 0 : 1;

  if(propError == 1){
    printf("ERROR! Move proposals should sum to 1.0\n");
    exit(EXIT_FAILURE);
  }

  double toleranceSq = tolerance*tolerance;
  double epsilon = EPSILON/AVOGADRO_CNST;
  double DIM_TEMPS[20], TEMPS[20];
  
  //DIM_TEMPS[0] = DIM_TEMP;
  //TEMPS[0] = DIM_TEMP*K_B/epsilon;

  /*USE REDUCED TEMPERATURES*/
  //TEMPS[0] = DIM_TEMP;  

  if((inputPtr=fopen(tempsFile,"r"))==NULL){
    printf("Cannot read temperatures from file %s\n", tempsFile);
    exit(EXIT_FAILURE);
  }

  int i;
  for(i=0;i<NUM_TEMPS;i++){
    //DIM_TEMPS[i] = RAT*DIM_TEMPS[i-1];
    //TEMPS[i] = RAT*TEMPS[i-1];
    fscanf(inputPtr,"%lf",&TEMPS[i]);
  }
  fclose(inputPtr);

  int doSwaps = 1;
  if(NUM_TEMPS == 1)
    doSwaps = 0;

  printf("N: %d\n", N);
  printf("LAMBDA: %f\n", LAMBDA);
  //printf("TEMPS (Kelvin): [%.2f", DIM_TEMPS[0]);
  //for(i=1;i<NUM_TEMPS;i++)
    //printf(", %.2f", DIM_TEMPS[i]);
  //printf("]\n");
  printf("TEMPS (Dimensionless): [%f",TEMPS[0]);
  for(i=1;i<NUM_TEMPS;i++)
    printf(", %.4f",TEMPS[i]);
  printf("]\n");
  printf("TEMP_INTEREST: %d\n", TEMP_INTEREST);
  printf("sigma: %f\n", sigma);
  printf("stretching region: [%d, %d]\n",stretchStart,stretchFinish);
  printf("PHI_MAX: %f\n", PHI_MAX);
  printf("fileInput: %d\n", fileInput);
  printf("configsInput: %s\n", configsInput);
  printf("moveSizesDir: %s\n",moveSizesDir);
  printf("Sampling window limits: [%d, %d]\n", Nmin, Nmax);
  printf("kappa: %f\n", kappa);
  printf("kappaRHS: %f\n", kappaRHS);
  //printf("updateInfoDir: %s\n", updateInfoDir);
  printf("update: %d\n", update);
  printf("tolerance: %f\n", tolerance);
  printf("PRINT_INTERVAL: %.2e\n", PRINT_INTERVAL);
  printf("PHI_UPDATE_INTERVAL: %.2e\n", PHI_UPDATE_INTERVAL);
  printf("UPDATE_INTERVAL: %.2e\n", UPDATE_INTERVAL*SWAP_INTERVAL);
  printf("SWAP_INTERVAL: %.2e\n", SWAP_INTERVAL);
  printf("NUM_SWAPS: %d\n", NUM_SWAPS);
  printf("NUM_SWAP_INTERVALS: %.2e\n", NUM_SWAP_INTERVALS);
  printf("RESET_INTERVAL: %.2e\n", RESET_INTERVAL);
  printf("SMALL_ANGLE_MOVES: %.2f\n", SMALL_ANGLE_MOVES);
  printf("REPTATION_MOVES: %.2f\n", REPTATION_MOVES);
  printf("CRANK_MOVES: %.2f\n", CRANK_MOVES);
  printf("END_ROT_MOVES: %.2f\n", END_ROT_MOVES);
  printf("END_BRIDGE_MOVES: %.2f\n", END_BRIDGE_MOVES);
  printf("Output Directory: %s\n", outputDir);
  printf("weightFile: %s\n", weightFile);
  printf("initCoords File: %s\n", initCoords);

  //int STATES = N+1;
  int STATES = 1000;
  double x_all[60][500], x_temp[20][500];
  double torsionAngles_all[60][500];
  double bondAngles_all[60][500];
  double bondAngles_temp[500]={0.0}, torsionAngles_temp[500]={0.0};
  long unsigned int occWeighted[20][1000]={0};
  //long unsigned int energyOccWeighted[20][1000];
  int j, k, swapCtr;
  long unsigned int seed[20];

  struct drand48_data randBuffer[20];
  for(i=0;i<NUM_TEMPS;i++){
    seed[i] = time(NULL) + i;
    //seed[i] = 100;
    srand48_r(seed[i],&randBuffer[i]);    
  }    

  struct drand48_data swapBuffer;
  long unsigned int swapSeed = time(NULL)+100;
  srand48_r(swapSeed,&swapBuffer);
  
  //long int seed = -100;
  double phiMax[20][1000], crankMax[20][1000],endMax[20];
  //  int particleState[N];
  int weightIndices[1000];
  double weightREAD[1000], weight[1000]={0.0};
  double weightSum[20]={0.0};
  FILE *weightPtr;
  double Eswap[20];
  double deltEswap;
  int swapped[20], allSwapped, firstSwapTemp, secondSwapTemp, swapCount[20]={0},swapCountA[20]={0}, doSwap;
  int accepted;  

  double occUnbiased[1000], freeEnergy[1000];  
  char occBiasedOutputFile[100], freeEnergyOutputFile[100];

  //for tracking zero crossings
  int rept[20];
  int reptZero[20], crossingCount[20];

  //for tracking move acceptance ratios
  long unsigned int accSmlAngleMoves[20]={0}, accSmlAngleMovesA[20]={0};
  long unsigned int accBigAngleMoves[20]={0}, accBigAngleMovesA[20]={0};
  long unsigned int accReptMoves[20]={0}, accReptMovesA[20]={0};
  long unsigned int accCrankMoves[20]={0}, accCrankMovesA[20]={0};
  long unsigned int accEndRotMoves[20]={0}, accEndRotMovesA[20]={0};
  long unsigned int accEndBrMoves[20]={0}, accEndBrMovesA[20]={0};


  int dummyInt;
  int config_tracking[20];

  //for logging configurations
  int printed[1000]={0}; 

  //for initialising with input files
  FILE *initConfigsPtr, *initMoveSizesPtr;
  char moveSizesFile[100];
  char restartInfoFile[100];
  FILE *restartInfoPtr;

  //for regular output of occupancies
  FILE *occOutputPtr;
  char occOutputFile[100];
  int resetCount = 0;
  int occWeightedSum;

/*  
  //for tracking swaps  
  char swapTrackingFile[100];  
  FILE *swapTrackingPtr;
  

  for(i=0;i<NUM_TEMPS;i++){
    config_tracking[i] = i;
    sprintf(swapTrackingFile,"%s/swapTracking/trace_chain%d.dat",outputDir,i);
    swapTrackingPtr = fopen(swapTrackingFile,"w");
    fprintf(swapTrackingPtr,"%d\n",i);
    fclose(swapTrackingPtr);
  }
  */

  initialise(x_all,N,initCoords);

  //start with a perfect hcp crystal
  //hcpGen(x_all,N);

  for(i=0;i<NUM_TEMPS;i++){
    for(j=0;j<3;j++){
      for(k=0;k<N;k++)
	x_all[j+3*i][k] = x_all[j][k];
    }
    endMax[i] = PHI_MAX;
  }
  weightPtr = fopen(weightFile,"r");
  //for(i=0;i<STATES;i++){
  int lines=0;
  while(!feof(weightPtr)){
    if(lines<STATES){
      fscanf(weightPtr,"%d %lf",&weightIndices[lines],&weightREAD[lines]);
      weightREAD[lines] *= TEMPS[TEMP_INTEREST];
      lines++;
    }
    else
      break;
  }
  fclose(weightPtr);
  for(i=0;i<lines-1;i++)
    weight[weightIndices[i]] = weightREAD[i];

  positions2bonds(x_all, bondAngles_temp, torsionAngles_temp, N);

  for(i=0;i<=N;i++){
    for(j=0;j<NUM_TEMPS;j++){
      bondAngles_all[j][i] = bondAngles_temp[i];
      torsionAngles_all[j][i] = torsionAngles_temp[i];
      occWeighted[j][i] = 0;
      phiMax[j][i] = PHI_MAX;
      crankMax[j][i] = PHI_MAX;
    }
  }

  if(fileInput==1){
    initConfigsPtr = fopen(configsInput,"r");
    for(i=0;i<N;i++){
      for(j=0;j<NUM_TEMPS;j++){
	for(k=0;k<3;k++)
	  fscanf(initConfigsPtr,"%lf",&x_all[3*j+k][i]);
      }
    }
    fclose(initConfigsPtr);
    sprintf(moveSizesFile,"%sphiMax.dat",moveSizesDir);
    initMoveSizesPtr = fopen(moveSizesFile,"r");
    for(i=0;i<N;i++){
      for(j=0;j<NUM_TEMPS;j++)
	fscanf(initMoveSizesPtr,"%lf",&phiMax[j][i]);
    }
    fclose(initMoveSizesPtr);

    sprintf(moveSizesFile,"%scrankMax.dat",moveSizesDir);
    initMoveSizesPtr = fopen(moveSizesFile,"r");
    for(i=0;i<N;i++){
      for(j=0;j<NUM_TEMPS;j++)
	fscanf(initMoveSizesPtr,"%lf",&crankMax[j][i]);
    }
    fclose(initMoveSizesPtr);

    sprintf(moveSizesFile,"%sendMax.dat",moveSizesDir);
    initMoveSizesPtr = fopen(moveSizesFile,"r");
    for(i=0;i<NUM_TEMPS;i++)
      fscanf(initMoveSizesPtr,"%lf",&endMax[i]);

    fclose(initMoveSizesPtr);
    

  }


  /*
  for(i=0;i<1000;i++){
    for(j=0;j<NUM_TEMPS;j++)
      energyOccWeighted[j][i]=0;
  }
  */
  for(i=0;i<STATES;i++){
    if(i<Nmin)
      weight[i] = -0.5*TEMPS[TEMP_INTEREST]*kappa*(1.0*(i-Nmin)*(i-Nmin))+weight[i];
    else{
      if(i>Nmax)
	weight[i] = -0.5*TEMPS[TEMP_INTEREST]*kappaRHS*(1.0*(i-Nmax)*(i-Nmax))+weight[i];
    }
  }

  FILE *weightOutputPtr;
  char weightOutputFile[100];
  sprintf(weightOutputFile,"%sweightFull.dat",outputDir);
  weightOutputPtr = fopen(weightOutputFile,"w");
  for(i=0;i<STATES;i++)
    fprintf(weightOutputPtr,"%d %f\n",i,weight[i]/TEMPS[TEMP_INTEREST]);
  fclose(weightOutputPtr);

  for(i=0;i<NUM_TEMPS;i++){
    //weightSum[i] = 0.0;
    rept[i] = 0;
    crossingCount[i] = 0;
    reptZero[i] = 0;
  }

  //for(i=0;i<=Nmax+1;i++)
  //printf("weight[%d]: %f\n", i, weight[i]);

  int SWAP_ATTEMPTS = 0;
  int firstParse[20];
  for(i=0;i<NUM_TEMPS;i++)
    firstParse[i] = 1;
  omp_set_num_threads(NUM_TEMPS);

  FILE *EtracePtr;
  char EtraceFile[100];
  int Etrace;
  Etrace = -E_swSum(x_all,N,LAMBDA);
  for(i=0;i<NUM_TEMPS;i++){
    sprintf(EtraceFile,"%sEtrace%d.dat",outputDir,i);
    EtracePtr = fopen(EtraceFile,"w");
    fprintf(EtracePtr,"%d\n", Etrace);
    fclose(EtracePtr);
  }

  FILE *updatePtr, *newBiasPtr, *kappaPtr;
  int fileUpdate;
  char updateYesNoFile[100], kappaFile[100], newBiasFile[100];
  //sprintf(updateYesNoFile,"%supdateYesNo.dat",updateInfoDir);
  //sprintf(kappaFile,"%snewKappa.dat",updateInfoDir);
  //sprintf(newBiasFile,"%snewBias.dat",updateInfoDir);

  //printf("updateYesNoFile: %s\n", updateYesNoFile);
  //printf("kappaFile: %s\n", kappaFile);
  //printf("newBiasFile: %s\n", newBiasFile);
 
  double bestEtotal[1000], worstEtotal[1000], bestE_LJcryst[1000], worstE_LJcryst[1000];
  double bestE_LJnonCryst[1000], worstE_LJnonCryst[1000];
  double bestTorsionPrint[1000], worstTorsionPrint[1000], bestE_totalPrint[1000], worstE_totalPrint[1000];
  double avEtotal[1000], avE_LJcryst[1000], avE_LJnonCryst[1000], avTorsionPrint[1000];
  double EtotalSum[1000], E_LJcrystSum[1000], E_LJnonCrystSum[1000], torsionPrintSum[1000];

  long int upwardMovesAttempted[20][1000]={0};
  long int upwardMovesAccepted[20][1000]={0};

  for(i=0;i<N;i++){
    bestEtotal[i] = 1000.0;
    bestE_LJcryst[i] = 1000.0;
    worstEtotal[i] = -1000.0;
    worstE_LJcryst[i] = -1000.0;
    EtotalSum[i] = 0.0;
    E_LJcrystSum[i] = 0.0;
    E_LJnonCrystSum[i] = 0.0;
    torsionPrintSum[i] = 0.0;        
  }
  
  double random_number;// random_sum=0.0;

  double run_time;
  run_time = omp_get_wtime();
  int equilibrated = 1;
  do{
#pragma omp parallel    									\
      default(none)							\
  shared(weightSum,upwardMovesAttempted,upwardMovesAccepted,swapCount,kappa,printed,randBuffer,seed,firstParse,rept,crossingCount,reptZero,update, occWeighted,Nmin,Nmax,PRINT_INTERVAL,PHI_UPDATE_INTERVAL,UPDATE_INTERVAL,ALPHA,BETA,GAMMA,sigma,stretchStart,stretchFinish,SMALL_ANGLE_MOVES,REPTATION_MOVES,CRANK_MOVES,END_ROT_MOVES,END_BRIDGE_MOVES,SWAP_INTERVAL,SWAP_ATTEMPTS,NUM_TEMPS,TEMPS,STATES,N,LAMBDA,LAMBDA_SQ,x_all,bondAngles_all,torsionAngles_all,phiMax,crankMax,endMax,weight,Eswap,outputDir,VMDdir,bestEtotal,worstEtotal,bestE_LJcryst,worstE_LJcryst,bestE_LJnonCryst,worstE_LJnonCryst,bestTorsionPrint,worstTorsionPrint,bestE_totalPrint,worstE_totalPrint,accSmlAngleMoves,accSmlAngleMovesA,accBigAngleMoves,accBigAngleMovesA,accReptMoves,accReptMovesA,accCrankMoves,accCrankMovesA,accEndRotMoves,accEndRotMovesA,accEndBrMoves,accEndBrMovesA, equilibrated,avEtotal,avE_LJcryst,avE_LJnonCryst,avTorsionPrint,EtotalSum,E_LJcrystSum,E_LJnonCrystSum,torsionPrintSum) \
  private(i,j,k,occBiasedOutputFile,freeEnergyOutputFile,random_number)	\
  firstprivate(occUnbiased,freeEnergy) 
    {
      int chain = omp_get_thread_num();
     
      double x[60][500];
      double x_old[60][500];
      double torsionAngles[500];
      double torsionAnglesNEW[500];
      double bondAngles[500];
      double bondAnglesNEW[500];
      int particleState[500];
      //double weightSum = 0.0;
      double beta = 1.0/TEMPS[chain];
      int cmpt;
      int nthBond, accepted;
      int move_type;
      //double E_LJ, E_LJpartial_BEF, E_LJpartial_AFT, E_torsion, newE_torsion;
      double E_sw, newE_sw, E_bond, newE_bond, E_torsion, newE_torsion;
      //double deltE_LJ, deltE_torsion, deltEtotal;
      double deltE_sw, deltE_bond, deltE_torsion, deltEtotal;
      double E_stretch, newE_stretch, deltE_stretch;
      double phiNew, thetaCrank, thetaNew; //phi_rept;
      //int acMovesTot;
      int N_crystalOld, N_crystalNew;
      int current_OrderParameter;
      long unsigned int occWeightedSum=1, occUniSum;
      double weightOld, weightNew, deltWeight, avWeight, OccUni, diff, accFactor;      
      //double energyFreeEnergy[1000];      
      FILE *occBiasedOutputPtr, *freeEnergyOutputPtr;
      //FILE *energyOccOutputPtr, *energyFreeEnergyOutputPtr;      
      FILE *crossingOutputPtr;
      char crossingOutputFile[100];

      int acMoves[500]={0}, acMovesA[500]={0}, acMovesCR[500]={0}, acMovesCRA[500]={0};
      int acMovesENDROT=0, acMovesENDROTA=0;
      double accRatio[500], accRatioCR[500];

      double tempE;

      char printDir[100], energyFile[100];
      FILE *energyPtr;

      double rndnum_1, rndnum_2;
      double rndnum_1sq, rndnum_2sq;

      int changed;

      int FEzero, firstTime, nextResolved=-1;
      double linGrad, linConst;
      
      FILE *configPtr;     
      FILE *moveAttemptsPtr;
      char moveAttemptsFile[100];    

      if(firstParse[chain] == 1){
	sprintf(crossingOutputFile,"%scrossings_chain%d.dat",outputDir,chain);
	crossingOutputPtr = fopen(crossingOutputFile,"w");
	fprintf(crossingOutputPtr," ");
	fclose(crossingOutputPtr);
	firstParse[chain] = 0;
      }         
      
      for(i=0;i<N;i++){
	for(cmpt=0;cmpt<3;cmpt++)
	  x[cmpt][i] = x_all[cmpt+3*chain][i];
	torsionAngles[i] = torsionAngles_all[chain][i];
  bondAngles[i] = bondAngles_all[chain][i];	
      }
         
      //E_LJ = E_LJsum(x,N);
      E_sw = E_swSum(x,N,LAMBDA);
      E_bond = Ebond(bondAngles,N);
      E_torsion = computeE_torsion(torsionAngles,N,0,0);
      E_stretch = Estretch(x,N,sigma,stretchStart,stretchFinish);

      //N_crystalOld = orderParam(x,N,0,particleState);
      //N_crystalOld = DETorderParam(x,N_Cstar,N,0,0.0,0,0,particleState);
      N_crystalOld = SNorderParam(x,N,LAMBDA,particleState,E_sw);      
      //sprintf(printDir,"%sbest/",VMDdir);
      //VMDprint(x,N,N_crystalOld,particleState,printDir);
    
      for(int step=1;step<=SWAP_INTERVAL;step++) {
	
	for(cmpt=0;cmpt<3;cmpt++){
	  for(i=0;i<N;i++)
	    x_old[cmpt][i] = x[cmpt][i];
	}

	
	//E_LJpartial_BEF = E_LJsum(x_old,N);
	weightOld = weight[N_crystalOld];

	drand48_r(&randBuffer[chain], &random_number);

  if(random_number < SMALL_ANGLE_MOVES)
    move_type = 0; //pivot move
  else{
    if(random_number < SMALL_ANGLE_MOVES + REPTATION_MOVES){
      if(random_number < SMALL_ANGLE_MOVES + REPTATION_MOVES/2.0)
        move_type = 1; //forward reptation move      
      else
        move_type = 2; //backward reptation move        
    }
    else{
      if(random_number < SMALL_ANGLE_MOVES + REPTATION_MOVES + CRANK_MOVES)
        move_type = 3; //crank-shaft move
      else{
        if(random_number < SMALL_ANGLE_MOVES + REPTATION_MOVES + CRANK_MOVES + END_ROT_MOVES){
          if(random_number < SMALL_ANGLE_MOVES + REPTATION_MOVES + CRANK_MOVES + END_ROT_MOVES/2.0)
            move_type = 4; //end "0" rotation                 
          else
            move_type = 5; //end "N" rotation                  
        }
        else{
          if(random_number < SMALL_ANGLE_MOVES + REPTATION_MOVES + CRANK_MOVES + END_ROT_MOVES + END_BRIDGE_MOVES/2.0)
            move_type = 6;
          else
            move_type = 7;
        }
      }
    }
  }  

  //printf("MADE IT (chain %d)\n",chain);

  switch(move_type){
    case 0: //pivot move
      accSmlAngleMoves[chain]++;
      drand48_r(&randBuffer[chain], &random_number);      
      nthBond = (int)(random_number*(N-3))+2;
      acMoves[nthBond]++;        
      for(i=0;i<N;i++)
        torsionAnglesNEW[i] = torsionAngles[i];
      drand48_r(&randBuffer[chain], &random_number);      
      phiNew = (2.0*random_number-1.0)*phiMax[chain][nthBond];
      torsionAnglesNEW[nthBond] += phiNew;
      rotateAboutU(x,N,nthBond,phiNew);
      //newE_sw = E_swSum(x,N,LAMBDA);
      newE_sw = E_sw + pivotMove_DE(x,x_old, N, nthBond, LAMBDA_SQ);
      newE_bond = E_bond;
      newE_torsion = computeE_torsion(torsionAnglesNEW,N,0,0);          
      newE_stretch = Estretch(x,N,sigma,stretchStart,stretchFinish);
      N_crystalNew = SNorderParam(x,N,LAMBDA,particleState,newE_sw);      
      weightNew = weight[N_crystalNew];
      deltE_sw = newE_sw-E_sw;
      deltE_torsion = newE_torsion-E_torsion;
      deltE_stretch = newE_stretch-E_stretch;
      deltWeight = weightNew-weightOld;
      deltEtotal = beta*(ALPHA*deltE_sw+BETA*deltE_torsion+GAMMA*deltE_bond+deltE_stretch-deltWeight);    
      
      break;
    case 1: //forward reptation move
      accReptMoves[chain]++;
      rndnum_1sq = 1.0, rndnum_2sq = 1.0;      
      while(rndnum_1sq + rndnum_2sq > 1.0){
        drand48_r(&randBuffer[chain],&rndnum_1);
        drand48_r(&randBuffer[chain],&rndnum_2);        
        rndnum_1 = 2.0*rndnum_1-1.0;
        rndnum_2 = 2.0*rndnum_2-1.0;
        rndnum_1sq = rndnum_1*rndnum_1;
        rndnum_2sq = rndnum_2*rndnum_2;
      }
      reptation(x,N,1,rndnum_1,rndnum_2);
      //newE_sw = E_swSum(x,N,LAMBDA);
      newE_sw = E_sw + forwardRep_DE(x,x_old, N, LAMBDA_SQ);
      positions2bonds(x,bondAnglesNEW,torsionAnglesNEW,N);
      newE_bond = Ebond(bondAnglesNEW,N);
      newE_torsion = computeE_torsion(torsionAnglesNEW,N,0,0);
      newE_stretch = Estretch(x,N,sigma,stretchStart,stretchFinish);
      N_crystalNew = SNorderParam(x,N,LAMBDA,particleState,newE_sw);
      weightNew = weight[N_crystalNew];
      deltE_sw = newE_sw-E_sw;
      deltE_bond = newE_bond-E_bond;
      deltE_torsion = newE_torsion-E_torsion;
      deltE_stretch = newE_stretch-E_stretch;
      deltWeight = weightNew-weightOld;
      deltEtotal = beta*(ALPHA*deltE_sw+BETA*deltE_torsion+GAMMA*deltE_bond+deltE_stretch-deltWeight);      

      break;
    case 2: //backward reptation move
      accReptMoves[chain]++;
      rndnum_1sq = 1.0, rndnum_2sq = 1.0;
      while(rndnum_1sq + rndnum_2sq > 1.0){
        drand48_r(&randBuffer[chain],&rndnum_1);
        drand48_r(&randBuffer[chain],&rndnum_2);      
        rndnum_1 = 2.0*rndnum_1-1.0;
        rndnum_2 = 2.0*rndnum_2-1.0;
        rndnum_1sq = rndnum_1*rndnum_1;
        rndnum_2sq = rndnum_2*rndnum_2;
      }
      reptation(x,N,0,rndnum_1,rndnum_2);
      //newE_sw = E_swSum(x,N,LAMBDA);
      newE_sw = E_sw  +   backwardRep_DE(x,x_old, N, LAMBDA_SQ);
      positions2bonds(x,bondAnglesNEW,torsionAnglesNEW,N);
      newE_bond = Ebond(bondAnglesNEW,N);
      newE_torsion = computeE_torsion(torsionAnglesNEW,N,0,0);
      newE_stretch = Estretch(x,N,sigma,stretchStart,stretchFinish);
      N_crystalNew = SNorderParam(x,N,LAMBDA,particleState,newE_sw);
      weightNew = weight[N_crystalNew];
      deltE_sw = newE_sw-E_sw;
      deltE_bond = newE_bond-E_bond;
      deltE_torsion = newE_torsion-E_torsion;
      deltE_stretch = newE_stretch-E_stretch;
      deltWeight = weightNew-weightOld;
      deltEtotal = beta*(ALPHA*deltE_sw+BETA*deltE_torsion+GAMMA*deltE_bond+deltE_stretch-deltWeight);      



      
      break;
    case 3: //crank-shaft move
      accCrankMoves[chain]++;
      drand48_r(&randBuffer[chain],&random_number);      
      nthBond = (int)(random_number*(N-3))+1.0;
      acMovesCR[nthBond]++;
      drand48_r(&randBuffer[chain],&random_number);      
      thetaCrank = (2.0*random_number-1.0)*crankMax[chain][nthBond];
      crankShaft(x,N,nthBond,thetaCrank);
      //newE_sw = E_swSum(x,N,LAMBDA);
      newE_sw = E_sw  + singleParticleMove_DE(x,x_old, N, nthBond, LAMBDA_SQ);
      positions2bonds(x,bondAnglesNEW,torsionAnglesNEW,N);
      newE_bond = Ebond(bondAnglesNEW,N);
      newE_torsion = computeE_torsion(torsionAnglesNEW,N,0,0);
      newE_stretch = Estretch(x,N,sigma,stretchStart,stretchFinish);
      N_crystalNew = SNorderParam(x,N,LAMBDA,particleState,newE_sw);
      weightNew = weight[N_crystalNew];
      deltE_sw = newE_sw-E_sw;
      deltE_bond = newE_bond-E_bond;
      deltE_torsion = newE_torsion-E_torsion;
      deltE_stretch = newE_stretch-E_stretch;
      deltWeight = weightNew-weightOld;      
      deltEtotal = beta*(ALPHA*deltE_sw+BETA*deltE_torsion+GAMMA*deltE_bond+deltE_stretch-deltWeight);

     
      break;
    case 4: //end "0" rotation
      accEndRotMoves[chain]++;
      acMovesENDROT++;
      rndnum_1sq = 1.0, rndnum_2sq = 1.0;      
      while(rndnum_1sq + rndnum_2sq > 1.0){
        drand48_r(&randBuffer[chain],&rndnum_1);
        drand48_r(&randBuffer[chain],&rndnum_2);        
        rndnum_1 = 2.0*rndnum_1-1.0;
        rndnum_2 = 2.0*rndnum_2-1.0;
        rndnum_1sq = rndnum_1*rndnum_1;
        rndnum_2sq = rndnum_2*rndnum_2;
      }
      drand48_r(&randBuffer[chain],&random_number);      
      thetaNew = (2.0*random_number-1.0)*endMax[chain];
      endRotate(x,N,0,rndnum_1,rndnum_2,thetaNew);
      //newE_sw = E_swSum(x,N,LAMBDA);
      newE_sw = E_sw  +  singleParticleMove_DE(x,x_old, N,0, LAMBDA_SQ);
      positions2bonds(x,bondAnglesNEW,torsionAnglesNEW,N);
      newE_bond = Ebond(bondAnglesNEW,N);
      newE_torsion = computeE_torsion(torsionAnglesNEW,N,0,0);
      newE_stretch = Estretch(x,N,sigma,stretchStart,stretchFinish);
      N_crystalNew = SNorderParam(x,N,LAMBDA,particleState,newE_sw);
      weightNew = weight[N_crystalNew];
      deltE_sw = newE_sw-E_sw;
      deltE_bond = newE_bond-E_bond;
      deltE_torsion = newE_torsion-E_torsion;
      deltE_stretch = newE_stretch-E_stretch;
      deltWeight = weightNew-weightOld;      
      deltEtotal = beta*(ALPHA*deltE_sw+BETA*deltE_torsion+GAMMA*deltE_bond+deltE_stretch-deltWeight);

      
      break;
    case 5: //end "N" rotation
      accEndRotMoves[chain]++;
      acMovesENDROT++;
      rndnum_1sq = 1.0, rndnum_2sq = 1.0;      
      while(rndnum_1sq + rndnum_2sq > 1.0){
        drand48_r(&randBuffer[chain],&rndnum_1);
        drand48_r(&randBuffer[chain],&rndnum_2);        
        rndnum_1 = 2.0*rndnum_1-1.0;
        rndnum_2 = 2.0*rndnum_2-1.0;
        rndnum_1sq = rndnum_1*rndnum_1;
        rndnum_2sq = rndnum_2*rndnum_2;
      }
      drand48_r(&randBuffer[chain],&random_number);      
      thetaNew = (2.0*random_number-1.0)*endMax[chain];
      endRotate(x,N,1,rndnum_1, rndnum_2,thetaNew);
      //newE_sw = E_swSum(x,N,LAMBDA);
      newE_sw = E_sw  +  singleParticleMove_DE(x,x_old, N,N-1, LAMBDA_SQ);
      positions2bonds(x,bondAnglesNEW,torsionAnglesNEW,N);
      newE_bond = Ebond(bondAnglesNEW,N);
      newE_torsion = computeE_torsion(torsionAnglesNEW,N,0,0);
      newE_stretch = Estretch(x,N,sigma,stretchStart,stretchFinish);
      N_crystalNew = SNorderParam(x,N,LAMBDA,particleState,newE_sw);
      weightNew = weight[N_crystalNew];
      deltE_sw = newE_sw-E_sw;
      deltE_bond = newE_bond-E_bond;
      deltE_torsion = newE_torsion-E_torsion;
      deltE_stretch = newE_stretch-E_stretch;
      deltWeight = weightNew-weightOld;      
      deltEtotal = beta*(ALPHA*deltE_sw+BETA*deltE_torsion+GAMMA*deltE_bond+deltE_stretch-deltWeight);
      
      
      break;
    case 6: //end "0" bridge
      accEndBrMoves[chain]++;
      drand48_r(&randBuffer[chain], &rndnum_1);
      drand48_r(&randBuffer[chain], &rndnum_2);
      accFactor  = endBridge(x,N,0,rndnum_1,rndnum_2,&changed );
      if(accFactor*accFactor > 0.0){
        //newE_sw = E_swSum(x,N,LAMBDA);
	newE_sw = E_sw  +  singleParticleMove_DE(x,x_old, N,changed, LAMBDA_SQ);
        positions2bonds(x,bondAnglesNEW,torsionAnglesNEW,N);
        newE_bond = Ebond(bondAnglesNEW,N);
        newE_torsion = computeE_torsion(torsionAnglesNEW,N,0,0);
	newE_stretch = Estretch(x,N,sigma,stretchStart,stretchFinish);
        N_crystalNew = SNorderParam(x,N,LAMBDA,particleState,newE_sw);
        weightNew = weight[N_crystalNew];
        deltE_sw = newE_sw-E_sw;
        deltE_bond = newE_bond-E_bond;
        deltE_torsion = newE_torsion-E_torsion;
	deltE_stretch = newE_stretch-E_stretch;
        deltWeight = weightNew-weightOld;         
	deltEtotal = beta*(ALPHA*deltE_sw+BETA*deltE_torsion+GAMMA*deltE_bond+deltE_stretch-deltWeight)-accFactor;
      }
      else deltEtotal = 100.0;
      break;
    case 7: //end "N" bridge
      accEndBrMoves[chain]++;
      drand48_r(&randBuffer[chain], &rndnum_1);
      drand48_r(&randBuffer[chain], &rndnum_2);
      accFactor  = endBridge(x,N,1,rndnum_1,rndnum_2, &changed);
      if(accFactor*accFactor > 0.0){
        //newE_sw = E_swSum(x,N,LAMBDA);
	newE_sw = E_sw  +  singleParticleMove_DE(x,x_old, N,changed, LAMBDA_SQ);
        positions2bonds(x,bondAnglesNEW,torsionAnglesNEW,N);
        newE_bond = Ebond(bondAnglesNEW,N);
        newE_torsion = computeE_torsion(torsionAnglesNEW,N,0,0);
	newE_stretch = Estretch(x,N,sigma,stretchStart,stretchFinish);
        N_crystalNew = SNorderParam(x,N,LAMBDA,particleState,newE_sw);
        weightNew = weight[N_crystalNew];
        deltE_sw = newE_sw-E_sw;
        deltE_bond = newE_bond-E_bond;
        deltE_torsion = newE_torsion-E_torsion;
	deltE_stretch = newE_stretch-E_stretch;
        deltWeight = weightNew-weightOld;                        
	deltEtotal = beta*(ALPHA*deltE_sw+BETA*deltE_torsion+GAMMA*deltE_bond+deltE_stretch-deltWeight)-accFactor;
      }
      else deltEtotal = 100.0;
      break;
    default: //invalid move
      printf("ERROR! Invalid move option\n");
      exit(EXIT_FAILURE);    
  }      
  
  if(N_crystalNew>N_crystalOld)
    upwardMovesAttempted[chain][N_crystalOld]++;

	accepted = 0;
	if(deltEtotal < 75.0){
	  if(deltEtotal <= 0.0)
	    accepted = 1;
	  else{
	    drand48_r(&randBuffer[chain], &random_number);            
	    if(exp(-deltEtotal) > random_number)
	      accepted = 1;
	  }
	}
	if(accepted == 1){
	  if(N_crystalNew>N_crystalOld)
	    upwardMovesAccepted[chain][N_crystalOld]++;
    //energyOccWeighted[chain][-(int)newE_sw]++;
    switch(move_type){      
      case 0:
        accSmlAngleMovesA[chain]++;
        acMovesA[nthBond]++;
        break;
      case 1:
        accReptMovesA[chain]++;
        rept[chain]++;
        break;
      case 2:
        accReptMovesA[chain]++;
        rept[chain]--;
        break;
      case 3:
        accCrankMovesA[chain]++;
        acMovesCRA[nthBond]++;
        break;
      case 4:
        accEndRotMovesA[chain]++;
        acMovesENDROTA++;
        break;
      case 5:
        accEndRotMovesA[chain]++;
        acMovesENDROTA++;
        break; 
      case 6:
        accEndBrMovesA[chain]++;
        break;
      case 7:
        accEndBrMovesA[chain]++;
        break;             
    }
    for(i=0;i<N;i++){
      bondAngles[i] = bondAnglesNEW[i];
      torsionAngles[i] = torsionAnglesNEW[i];
    }
    E_sw = newE_sw;
    E_bond = newE_bond;
    E_torsion = newE_torsion;
    E_stretch = newE_stretch;
    N_crystalOld = N_crystalNew;
    weightOld = weightNew;
  }
  else{
    //energyOccWeighted[chain][-(int)E_sw]++;
    for(cmpt=0;cmpt<3;cmpt++){
      for(i=0;i<N;i++)
        x[cmpt][i] = x_old[cmpt][i];
    }
  }
	
	
  if(accepted ==1 && newE_sw>0.0){
    printf("move_type: %d\n", move_type);
    printf("accepted: %d\n", accepted);
    printf("N_crystalOld: %d -> N_crystalNew: %d\n", N_crystalOld,N_crystalNew); 
    printf("weightOld: %f -> weightNew: %f\n", weightOld,weightNew);
    printf("deltE_total: %f\n", deltEtotal); 
    printf("E_sw: %f -> newE_sw: %f\n", E_sw,newE_sw);
    printf("step: %d\n", step);
    printf("-----------------------------\n");
    exit(EXIT_FAILURE);
  }

  weightSum[chain] += 1.0/exp(beta*weightOld);
  occWeighted[chain][N_crystalOld]++;     

	if(bondLengthCheck(x,N)==1){
	  printf("***Bond Length Error!***\n");
	  printf("chain: %d\n", chain);
	  printf("nthBond: %d\n", nthBond);
	  printf("move_type: %d\n",move_type);
	  printf("step: %d\n", step);
	  printf("Input Coordinates\n");
	  printf("--------------------\n");
	  for(i=0;i<N;i++)
	    printf("(%f, %f, %f)\n", x[0][i],x[1][i],x[2][i]);
	  exit(EXIT_FAILURE);
	}
		      
	if(rept[chain]%N==0 && rept[chain]!=reptZero[chain]){
	  #pragma omp critical
	  {
	  crossingCount[chain]++;
	  reptZero[chain] = rept[chain];
	  crossingOutputPtr = fopen(crossingOutputFile,"a");
	  fprintf(crossingOutputPtr,"Crossing %d after %d steps\n",crossingCount[chain],(int)(SWAP_ATTEMPTS*SWAP_INTERVAL)+step);
	  fclose(crossingOutputPtr);
	  }
	}	
	
	if(chain == 0 ){
#pragma omp critical
	  {	             
	    //current_OrderParameter = orderParam(x,N,0,particleState);
      //current_OrderParameter = DETorderParam(x,N_Cstar,N,0,0.0,0,0,particleState);
      current_OrderParameter = SNorderParam(x,N,LAMBDA,particleState,E_sw);
	    if(N_crystalOld != current_OrderParameter){
        SNorderParam(x,N,LAMBDA,particleState,E_sw);
        sprintf(printDir,"%sbest/",VMDdir);        
        printf("move_type: %d\n", move_type);
        printf("accepted: %d\n", accepted);
        printf("current_OrderParameter: %d\n", current_OrderParameter);
        printf("N_crystalOld: %d\n", N_crystalOld);        
        printf("N_crystalNew: %d\n", N_crystalNew);        
		    printf("ERROR: Order Parameter Inconsistent\n");
		    exit( EXIT_FAILURE);
	    }
	    
	    /*
      if(printed[N_crystalOld] < 1 && accepted==1){
        printed[N_crystalOld]++;
        sprintf(printDir,"%scoords/state%d_%d.dat",VMDdir,N_crystalOld,printed[N_crystalOld]);
        configPtr = fopen(printDir,"w");
        for(i=0;i<N;i++)
          fprintf(configPtr,"%f %f %f\n", x[0][i],x[1][i],x[2][i]);
        fclose(configPtr);
        sprintf(printDir,"%svis/",VMDdir);
        VMDprint(x,N,N_crystalOld,printed[N_crystalOld],particleState,printDir);
      }
	    */

	    /*
	    tempE = ALPHA*E_sw+BETA*E_torsion+GAMMA*E_bond;	   

	    EtotalSum[N_crystalOld] += tempE;
	    avEtotal[N_crystalOld] = EtotalSum[N_crystalOld]/(1.0*occWeighted[0][N_crystalOld]);
	    //E_LJcrystSum[N_crystalOld] += ALPHA*E_LJcrystAv(x,N,particleState);
	    //avE_LJcryst[N_crystalOld] = E_LJcrystSum[N_crystalOld]/(1.0*occWeighted[0][N_crystalOld]);
	    //E_LJnonCrystSum[N_crystalOld] += ALPHA*E_LJnonCrystAv(x,N,particleState);
	    //avE_LJnonCryst[N_crystalOld] = E_LJnonCrystSum[N_crystalOld]/(1.0*occWeighted[0][N_crystalOld]);
	    torsionPrintSum[N_crystalOld] += BETA*E_torsion;
	    avTorsionPrint[N_crystalOld] = torsionPrintSum[N_crystalOld]/(1.0*occWeighted[0][N_crystalOld]);

	    if(tempE<bestEtotal[N_crystalOld]){
	      bestEtotal[N_crystalOld] = tempE;
	      //bestE_LJcryst[N_crystalOld] = ALPHA*E_LJcrystAv(x,N,particleState);
        //bestE_LJnonCryst[N_crystalOld] = ALPHA*E_LJnonCrystAv(x,N,particleState);
	      bestTorsionPrint[N_crystalOld] = BETA*E_torsion;
	      bestE_totalPrint[N_crystalOld] = ALPHA*E_sw+BETA*E_torsion+GAMMA*E_bond;
	      //sprintf(printDir,"%sbest/",VMDdir);
	      //VMDprint(x,N,N_crystalOld,particleState,printDir);
	      sprintf(energyFile,"%sbestEnergies.dat",printDir);
	      energyPtr = fopen(energyFile,"w");
	      fprintf(energyPtr,"n  E_LJcryst/part E_LJnoncryst/part   E_torsion      E_total\n");
	      fprintf(energyPtr,"------------------------------------------------------------\n");
	      for(i=0;i<N;i++){
		if(bestEtotal[i]<1000.0)
      if(i == 0)
        fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, 0.0, 0.0, bestTorsionPrint[i],bestE_totalPrint[i]);
        //fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, bestE_LJcryst[i],bestE_LJnonCryst[i]/(1.0*N), bestTorsionPrint[i],bestE_totalPrint[i]);
      else
        fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, 0.0, 0.0, bestTorsionPrint[i],bestE_totalPrint[i]);
		    //fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, bestE_LJcryst[i]/(1.0*i),bestE_LJnonCryst[i]/(1.0*(N-i)), bestTorsionPrint[i],bestE_totalPrint[i]);
	      }
	      fclose(energyPtr);
	    }
	    else{
	      if(tempE>worstEtotal[N_crystalOld]){
		worstEtotal[N_crystalOld] = tempE;
		//worstE_LJcryst[N_crystalOld] = ALPHA*E_LJcrystAv(x,N,particleState);
    //worstE_LJnonCryst[N_crystalOld] = ALPHA*E_LJnonCrystAv(x,N,particleState);
		worstTorsionPrint[N_crystalOld] = BETA*E_torsion;
		worstE_totalPrint[N_crystalOld] = ALPHA*E_sw+BETA*E_torsion+GAMMA*E_bond;
		sprintf(printDir,"%sworst/",VMDdir);
		//VMDprint(x,N,N_crystalOld,particleState,printDir);
		//sprintf(energyFile,"%sworstEnergies.dat",printDir);
		energyPtr = fopen(energyFile,"w");
		fprintf(energyPtr,"n  E_LJcryst/part E_LJnoncryst/part   E_torsion      E_total\n");
        fprintf(energyPtr,"------------------------------------------------------------\n");
		for(i=0;i<N;i++){
		  if(worstEtotal[i]>-1000.0)
		    if(i == 0)
          fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, 0.0, 0.0, bestTorsionPrint[i],bestE_totalPrint[i]);
        //fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, bestE_LJcryst[i],bestE_LJnonCryst[i]/(1.0*N), bestTorsionPrint[i],bestE_totalPrint[i]);
      else
        fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, 0.0, 0.0, bestTorsionPrint[i],bestE_totalPrint[i]);
        //fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, bestE_LJcryst[i]/(1.0*i),bestE_LJnonCryst[i]/(1.0*(N-i)), bestTorsionPrint[i],bestE_totalPrint[i]);
    }
		fclose(energyPtr);
	      }
	    }
      */
	  }	  
	}	

	if(step%(int)PHI_UPDATE_INTERVAL==0){
	  for(i=1;i<N;i++){
	    accRatio[i] = (1.0*acMovesA[i])/(1.0*acMoves[i]);
      accRatioCR[i] = (1.0*acMovesCRA[i])/(1.0*acMovesCR[i]);
	    acMoves[i] = 0;
	    acMovesA[i] = 0;
      acMovesCR[i] = 0;
      acMovesCRA[i] = 0;            
	    if(accRatio[i] < 0.5)
	      phiMax[chain][i] *= 0.95;
	    else
	      if( phiMax[chain][i]< PI/3.0)
		      phiMax[chain][i] *= 1.05;
      if(accRatioCR[i] < 0.5)
        crankMax[chain][i] *= 0.95;
      else
        if(crankMax[chain][i] < PI/3.0)
          crankMax[chain][i] *= 1.05;              
	  }
    
    if(((1.0*acMovesENDROTA)/(1.0*acMovesENDROT)) < 0.5){
      //if(endMax[chain] > 0.05)
        endMax[chain] *= 0.95;
    }
    else{
      if(endMax[chain] < PI/3.0)
        endMax[chain] *= 1.05;
    }    
    acMovesENDROT=0;
    acMovesENDROTA=0;
	}

	if(step%(int)PRINT_INTERVAL==0 && equilibrated == 1){
    /*
		if(chain == 0){    
		sprintf(printDir,"%saverage/",VMDdir);
		sprintf(energyFile,"%savEnergies.dat",printDir);
	    energyPtr = fopen(energyFile,"w");
	    fprintf(energyPtr,"n  E_LJcryst/part E_LJnoncryst/part   E_torsion      E_total\n");
	    fprintf(energyPtr,"------------------------------------------------------------\n");
	    for(i=0;i<N;i++){
	    if(occWeighted[0][i]>0){			
      		if(i == 0)
            fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, 0.0, 0.0, avTorsionPrint[i],avEtotal[i]);
        		//fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, avE_LJcryst[i],avE_LJnonCryst[i]/(1.0*N), avTorsionPrint[i],avEtotal[i]);
      		else
            fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, 0.0, 0.0, avTorsionPrint[i],avEtotal[i]);
		    	//fprintf(energyPtr,"%d %12f %17f %14f %14f\n", i, avE_LJcryst[i]/(1.0*i),avE_LJnonCryst[i]/(1.0*(N-i)), avTorsionPrint[i],avEtotal[i]);
	      }
	  }
	      fclose(energyPtr);
	  }
    */
		occWeightedSum=0;
	  for( i=0;i<STATES;i++)
	    occWeightedSum += occWeighted[chain][i];
	  avWeight = weightSum[chain]/(1.0*(occWeightedSum));	  
	  sprintf(freeEnergyOutputFile,"%sFE_chain%d.dat",outputDir,chain);
	  freeEnergyOutputPtr = fopen(freeEnergyOutputFile,"w");        
	  firstTime=0;
    nextResolved=-1;
	  for(i=0;i<=STATES;i++){	    	
	      occUnbiased[i] = ((1.0*occWeighted[chain][i])/exp(beta*weight[i]))/avWeight;
        if(occWeighted[chain][i] > THRESHOLD){	  
	         if(firstTime==0){
	           FEzero=i;
	           firstTime=1;
	         }	  
	         //freeEnergy[i] = log(occUnbiased[FEzero]/occUnbiased[i]);
		 freeEnergy[i] = -log(occUnbiased[i]);
        }
	      else{
          if(i==0)
            freeEnergy[Nmin] = 0.0;
          else{
            if(firstTime==1){
              for(j=i+1;j<STATES;j++){
                if(occWeighted[chain][j]>THRESHOLD){
                  nextResolved=j;
                  break;
                }
              }
              if(nextResolved>0){
                freeEnergy[nextResolved] = log(occUnbiased[FEzero]/occUnbiased[nextResolved]);
                linGrad = (freeEnergy[nextResolved]-freeEnergy[i-1])/(1.0*(nextResolved-i+1));
                linConst = freeEnergy[i-1]-linGrad*(i-1);
                for(j=i;j<nextResolved;j++)
                  freeEnergy[j] = linGrad*j+linConst;
                i=nextResolved;
                nextResolved=-1;
              }
              else
                freeEnergy[i] = freeEnergy[i-1];
            }
            else
              freeEnergy[i] = freeEnergy[i-1];
          }	       
	      }                
	
	/*
	    }	    
      else{
	if(occWeighted[chain][i] > THRESHOLD)
		      freeEnergy[i] = log((1.0*occWeighted[chain][Nmin])/(1.0*occWeighted[chain][i]))+weight[i]-weight[Nmin];
	      else{
		if(i>Nmin)
		  freeEnergy[i] = freeEnergy[i-1];
		else
		  freeEnergy[i] = 0.0;
	      }        
	    } 
	    */           

	if(occWeighted[chain][i] > THRESHOLD)
	      fprintf(freeEnergyOutputPtr,"%d %f\n", i, freeEnergy[i]);          
	  }
    /*
    sprintf(occBiasedOutputFile,"%senergyOccW_chain%d.dat",outputDir,chain);
    energyOccOutputPtr = fopen(occBiasedOutputFile,"w");
    for(i=0;i<1000;i++){
      if(energyOccWeighted[chain][i] > THRESHOLD){          
        energyFreeEnergy[i] = -log((1.0*energyOccWeighted[chain][i])/(1.0*occWeightedSum));  
        fprintf(energyFreeEnergyOutputPtr,"%d %f\n", -i,energyFreeEnergy[i]);
      }
      if(energyOccWeighted[chain][i] > 0)
        fprintf(energyOccOutputPtr,"%d %f\n", -i, (1.0*energyOccWeighted[chain][i])/(1.0*occWeightedSum));
    }
    */
	  fclose(freeEnergyOutputPtr);
    //fclose(energyFreeEnergyOutputPtr);
    //fclose(energyOccOutputPtr);
	  sprintf(occBiasedOutputFile,"%soccW_chain%d.dat",outputDir,chain);
	  occBiasedOutputPtr = fopen(occBiasedOutputFile,"w");

	  for(i=0;i<STATES;i++){
	    if(occWeighted[chain][i] > 0)
	      fprintf(occBiasedOutputPtr,"%d %f\n", i, (1.0*occWeighted[chain][i])/(1.0*(occWeightedSum)));
	  }    
	  fclose(occBiasedOutputPtr);

    //for(occUniSum=0, i=Nmin;i<=Nmax;i++)	   
	   // occUniSum += occWeighted[chain][i];

	  sprintf(moveAttemptsFile,"%supwardMoves_chain%d.dat",outputDir,chain);
	  moveAttemptsPtr = fopen(moveAttemptsFile,"w");	  
	  fprintf(moveAttemptsPtr,"State Up Moves Attempted Up Moves Accepted\n");
	  fprintf(moveAttemptsPtr,"------------------------------------------\n");
	  for(i=0;i<STATES;i++){
	    if(occWeighted[chain][i]>0){
	      fprintf(moveAttemptsPtr,"%3d %10f %10f\n", i,(1.0*upwardMovesAttempted[chain][i])/(1.0*occWeighted[chain][i]),(1.0*upwardMovesAccepted[chain][i])/(1.0*upwardMovesAttempted[chain][i]));
	    }
	  }
	  fclose(moveAttemptsPtr);
	}	
      }
      for(i=0;i<N;i++){
	for(cmpt=0;cmpt<3;cmpt++)
	  x_all[cmpt+3*chain][i] = x[cmpt][i];
  bondAngles_all[chain][i] = bondAngles[i];
	torsionAngles_all[chain][i] = torsionAngles[i];
      }
      Eswap[chain] = ALPHA*E_sw+BETA*E_torsion+GAMMA*E_bond-weight[N_crystalOld];

      #pragma omp barrier      
      //printf("Goodbye from chain %d %f %f %f \n",chain,x[0][0], x[0][N-1],torsionAnglesNEW[3] );  
    if(chain==0  && equilibrated == 1 && update == 1){
      if(occWeightedSum%(int)(UPDATE_INTERVAL*SWAP_INTERVAL)==0){	
	//if(occWeighted[chain][Nmin] > (int)1e4){
	  for(i=Nmin;i<=Nmax;i++)
	    weight[i] = TEMPS[0]*freeEnergy[i];
	  //}	  
	for(i=0;i<STATES;i++){
	  for(j=0;j<NUM_TEMPS;j++)
	    occWeighted[j][i] = 0;
	  EtotalSum[i] = 0.0;
   	  E_LJcrystSum[i] = 0.0;
   	  E_LJnonCrystSum[i] = 0.0;
   	  torsionPrintSum[i] = 0.0;	 
	}
	for(i=0;i<NUM_TEMPS;i++){
	  accSmlAngleMoves[i]=0;
	  accSmlAngleMovesA[i]=0;
	  accReptMoves[i]=0;
	  accReptMovesA[i]=0;
	  accEndRotMoves[i]=0;
	  accEndRotMovesA[i]=0;
	  accEndBrMoves[i]=0;
	  accEndBrMovesA[i]=0;
	}
  /*
  for(i=0;i<1000;i++){
    for(j=0;j<NUM_TEMPS;j++)
      energyOccWeighted[j][i] = 0;
  }
  */
	for(i=0;i<NUM_TEMPS;i++)
	  weightSum[i] = 0.0;
	UPDATE_INTERVAL = 1.0*floor(1.1*UPDATE_INTERVAL);
	/*
	kappaPtr = fopen(kappaFile,"r");
        fscanf(kappaPtr,"%lf",&kappa);
	fclose(kappaPtr);
	for(i=0;i<Nmin;i++)
	  weight[i] = -0.5*kappa*(1.0*(i-Nmin)*(i-Nmin));
	*/
	printf("WEIGHT FUNCTION UPDATED\n");
	printf("New update interval: %.2e\n", UPDATE_INTERVAL*SWAP_INTERVAL);
	//printf("kappa: %lf\n", kappa);
      }
      else{
        printf("*********************************\n");
        printf("* %7ld MC moves until update *\n",(int)(UPDATE_INTERVAL*SWAP_INTERVAL)-occWeightedSum);         
        printf("*********************************\n");
      }
    }
    } /*-----END OF PARALLEL REGION-----*/

    	  
    sprintf(restartInfoFile,"%srestartInfo/configs.dat",outputDir);
    restartInfoPtr = fopen(restartInfoFile,"w");
    for(i=0;i<N;i++){
      for(j=0;j<3*NUM_TEMPS;j++)
	fprintf(restartInfoPtr,"%f ",x_all[j][i]);
      fprintf(restartInfoPtr,"\n");
    }
    fclose(restartInfoPtr);
    sprintf(restartInfoFile,"%srestartInfo/phiMax.dat",outputDir);
    restartInfoPtr = fopen(restartInfoFile,"w");
    for(i=0;i<N;i++){
      for(j=0;j<NUM_TEMPS;j++)
	fprintf(restartInfoPtr,"%f ",phiMax[j][i]);
      fprintf(restartInfoPtr,"\n");
    }
    fclose(restartInfoPtr);
    sprintf(restartInfoFile,"%srestartInfo/crankMax.dat",outputDir);
    restartInfoPtr = fopen(restartInfoFile,"w");
    for(i=0;i<N;i++){
      for(j=0;j<NUM_TEMPS;j++)
	fprintf(restartInfoPtr,"%f ",crankMax[j][i]);
      fprintf(restartInfoPtr,"\n");
    }
    fclose(restartInfoPtr);
    sprintf(restartInfoFile,"%srestartInfo/endMax.dat",outputDir);
    restartInfoPtr = fopen(restartInfoFile,"w");
    for(i=0;i<NUM_TEMPS;i++)
      fprintf(restartInfoPtr,"%f\n",endMax[i]);
    fclose(restartInfoPtr);
    
    occWeightedSum=0;
    for(i=0;i<STATES;i++)
      occWeightedSum += occWeighted[0][i];
    if(occWeightedSum == RESET_INTERVAL){
      resetCount++;
      sprintf(occOutputFile,"%soccupancies/block_%d.dat",outputDir,resetCount);
      occOutputPtr = fopen(occOutputFile,"w");
      for(i=0;i<STATES;i++){
	for(j=0;j<NUM_TEMPS;j++)
	  fprintf(occOutputPtr,"%ld ",occWeighted[j][i]);
	fprintf(occOutputPtr,"\n");
      }
      fclose(occOutputPtr);
      for(i=0;i<STATES;i++){
      	for(j=0;j<NUM_TEMPS;j++)
      		occWeighted[j][i] = 0;
      }
      for (i=0;i<NUM_TEMPS;i++)
      	weightSum[i]=0.0;
      {
      	/* code */
      }
    }
    

    /******SERIAL CODE FOR ATTEMPTING SWAPS******/
  if(doSwaps == 1){ 
    for(k=0;k<NUM_TEMPS;k++){
      for(i=0;i<N;i++)
	for(j=0;j<3;j++)
	  x_temp[j][i] = x_all[j+3*k][i];
      Etrace = -E_swSum(x_temp,N,LAMBDA);
      sprintf(EtraceFile,"%sEtrace%d.dat",outputDir,k);      
      EtracePtr = fopen(EtraceFile,"a");
      fprintf(EtracePtr,"%d\n",Etrace);
      fclose(EtracePtr);
    }
    
    for(swapCtr=0;swapCtr<NUM_SWAPS;swapCtr++){
      drand48_r(&swapBuffer, &random_number);      
      firstSwapTemp = (int)((NUM_TEMPS-1)*random_number);
      secondSwapTemp = firstSwapTemp+1;
      swapCount[firstSwapTemp]++;      
      deltEswap = (1.0/TEMPS[firstSwapTemp]-1.0/TEMPS[secondSwapTemp])*(Eswap[secondSwapTemp]-Eswap[firstSwapTemp]);      
      accepted = 0;
      if(deltEswap < 75.0){
	if(deltEswap <= 0.0)
	  accepted = 1;
	else{
	  drand48_r(&swapBuffer, &random_number);           
	  if(exp(-deltEswap) > random_number)
	    accepted = 1;
	}
      }
      if(accepted == 1){	
	swapCountA[firstSwapTemp]++;
	for(i=0;i<N;i++){
	  for(j=0;j<3;j++){
	    x_temp[j][i] = x_all[j+3*firstSwapTemp][i];
	    x_all[j+3*firstSwapTemp][i] = x_all[j+3*secondSwapTemp][i];
	    x_all[j+3*secondSwapTemp][i] = x_temp[j][i];
	  }
	  torsionAngles_temp[i] = torsionAngles_all[firstSwapTemp][i];
	  torsionAngles_all[firstSwapTemp][i] = torsionAngles_all[secondSwapTemp][i];
	  torsionAngles_all[secondSwapTemp][i] = torsionAngles_temp[i];
	}
  
	dummyInt =  config_tracking[firstSwapTemp];
	config_tracking[firstSwapTemp] = config_tracking[secondSwapTemp];
	config_tracking[secondSwapTemp] = dummyInt;	  
      }
      for(k=0;k<NUM_TEMPS;k++){
	if(firstSwapTemp == k || secondSwapTemp == k){
	  for(i=0;i<N;i++)
	    for(j=0;j<3;j++)
	      x_temp[j][i] = x_all[j+3*k][i];
	  Etrace = -E_swSum(x_temp,N,LAMBDA);
	  sprintf(EtraceFile,"%sEtrace%d.dat",outputDir,k);      
	  EtracePtr = fopen(EtraceFile,"a");
	  fprintf(EtracePtr,"%d\n",Etrace);
	  fclose(EtracePtr);
	}
      }
    }
    
    /*
    if(SWAP_ATTEMPTS<1000){
    for(i=0;i<NUM_TEMPS;i++){
      sprintf(swapTrackingFile,"%s/swapTracking/trace_chain%d.dat",outputDir,i);
      swapTrackingPtr = fopen(swapTrackingFile,"a");
      fprintf(swapTrackingPtr,"%d\n", config_tracking[i]);
      fclose(swapTrackingPtr);
    }
	}
    */
 
    printf("Swap Acceptance Ratios: \n");
 
    int swapTot=0;
    //int swapTotA=0;
    for(i=0;i<NUM_TEMPS-1;i++){
      printf("Swaps from Temp %d to Temp %d: %f (%d out of %d)\n",i,i+1,(1.0*swapCountA[i])/(1.0*swapCount[i]),swapCountA[i],swapCount[i]);
      swapTot += swapCount[i];
      //swapTotA += swapCountA[i];
    }
    printf("Swap Tot: %d\n", swapTot);
    //printf("Swap Tot Acc: %d\n", swapTotA);
  }
  printf("Move Acceptance Ratios\n");
  printf("------------------------------\n");
  printf("TEMP        SML ANGLE   REPT          CRANK      END ROT          END BRIDGE\n");
  for(i=0;i<NUM_TEMPS;i++)
    //printf("%2d (%3dK) %10f %10f %10f %10f %10f\n", i, (int)DIM_TEMPS[i], (1.0*accSmlAngleMovesA[i])/(1.0*accSmlAngleMoves[i]), (1.0*accReptMovesA[i])/(1.0*accReptMoves[i]), (1.0*accCrankMovesA[i])/(1.0*accCrankMoves[i]),(1.0*accEndRotMovesA[i])/(1.0*accEndRotMoves[i]),(1.0*accEndBrMovesA[i])/(1.0*accEndBrMoves[i]));
    printf("%2d (%5.2f) %10f %10f %10f %10f %10f\n", i, TEMPS[i], (1.0*accSmlAngleMovesA[i])/(1.0*accSmlAngleMoves[i]), (1.0*accReptMovesA[i])/(1.0*accReptMoves[i]), (1.0*accCrankMovesA[i])/(1.0*accCrankMoves[i]),(1.0*accEndRotMovesA[i])/(1.0*accEndRotMoves[i]),(1.0*accEndBrMovesA[i])/(1.0*accEndBrMoves[i]));

  printf("Small Angle Move Sizes (Bond 3 and Bond N/2)  Crank-Shaft Size  End Rot\n");
  for(i=0;i<NUM_TEMPS;i++)
    printf("%2d (%5.2f) %10f %10f %10f %10f %10f\n", i, TEMPS[i], phiMax[i][3], phiMax[i][N/2], crankMax[i][3],crankMax[i][N/2],endMax[i]);

  long int occWeightedSum[20]={0};
  for( i=0;i<=STATES;i++){
    for(j=0;j<NUM_TEMPS;j++)
      occWeightedSum[j] += occWeighted[j][i];
  }

  printf("pivot moves: %f\n", 1.0*accSmlAngleMoves[0]/(1.0*occWeightedSum[0]));
  printf("rept moves: %f\n", 1.0*accReptMoves[0]/(1.0*occWeightedSum[0]));
  printf("crank-shaft moves: %f\n", 1.0*accCrankMoves[0]/(1.0*occWeightedSum[0]));
  printf("end moves: %f\n", 1.0*accEndRotMoves[0]/(1.0*occWeightedSum[0]));
  printf("end bridge moves: %f\n", 1.0*accEndBrMoves[0]/(1.0*occWeightedSum[0]));  
  
  equilibrated = 1;
 if(equilibrated == 0){
   	if(SWAP_ATTEMPTS<4){
   		printf("%d swap cycles until equilibrated\n",4-SWAP_ATTEMPTS);
   	}
   	else{
	  
	  //if(occWeighted[0][Nmin] > 1000){   		
   		  for(i=0;i<=N;i++){
   		   	for(j=0;j<NUM_TEMPS;j++)
   			  	occWeighted[j][i] = 0;
   			  EtotalSum[i] = 0.0;
   			  E_LJcrystSum[i] = 0.0;
   			  E_LJnonCrystSum[i] = 0.0;
   			  torsionPrintSum[i] = 0.0;
   		 }
   		 printf("***EQUILIBRATED***\n");
   		 equilibrated = 1;    		
		 // }
	//else
        //printf("Waiting for sufficient visits to Nmin\n");
   }
 }
  
  SWAP_ATTEMPTS++;
  
  /*
    for(i=0;i<5;i++)
    printf("(%f,%f,%f)\n", x_all[0][i], x_all[1][i], x_all[2][i]);
  */

  /*
  updatePtr = fopen(updateYesNoFile,"r");
  fscanf(updatePtr,"%d",&fileUpdate);
  fclose(updatePtr);
  if(fileUpdate == 1){
    kappaPtr = fopen(kappaFile,"r");
    fscanf(kappaPtr,"%lf",&kappa);
    fclose(kappaPtr);
    newBiasPtr = fopen(newBiasFile,"r");
    lines=0;
    while(!feof(newBiasPtr)){
      if(lines<STATES){
	fscanf(newBiasPtr,"%d %lf",&weightIndices[lines],&weightREAD[lines]);
	weightREAD[lines] *= TEMPS[TEMP_INTEREST];
	lines++;
      }
      else
	break;
    }
    fclose(newBiasPtr);
    for(i=0;i<lines-1;i++)
      weight[weightIndices[i]] = weightREAD[i];
    for(i=0;i<STATES;i++){
      if(i<Nmin)
	weight[i] = -0.5*kappa*(1.0*(i-Nmin)*(i-Nmin))+weight[Nmin];
      else
	if(i>Nmax)
	  weight[i] = -0.5*kappaRHS*(1.0*(i-Nmax)*(i-Nmax))+weight[Nmax];
    }
    for(i=0;i<NUM_TEMPS;i++){
      for(j=0;j<STATES;j++)
	occWeighted[i][j] = 0;
      weightSum[i] = 0.0;
    }
    printf("WEIGHT FUNCTION UPDATED FROM FILE: %s\n", newBiasFile);
  }
  */

  
  } while(SWAP_ATTEMPTS < NUM_SWAP_INTERVALS);

  for(i = 0  ;   i <= NUM_TEMPS-1    ;   i++){
    printf("Chain %d %f %f %f\n",i,x_all[3*i][0], x_all[3*i][N-1],torsionAngles_all[i][3] );
  }
  run_time = omp_get_wtime()-run_time;
  int hrs = (int)run_time / 3600;
  int mins = ((int)run_time-hrs*3600) / 60;
  double secs = run_time-hrs*3600-mins*60;
  printf("run_time: %d hours, %d mins, %f seconds\n", hrs, mins, secs);

      
  return 0;
}


