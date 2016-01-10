/*VMDprint.c
Script to create a VMD file with given coordinates to produce a 3D image of the system, with bonds.
*/

#include <stdio.h>
#include <stdlib.h>

#define SIZE 500

void VMDprint(double x[][500], int N, int k, int version,int *particleState,char *VMDdir)
{  

  FILE *coordsPtr, *bondsPtr, *vmdPtr;
  int totalBonds = N-1, i, j, bondCount = 0;
  double scaleFactor = 1.0;///0.382;
  int bondsMatrix[SIZE][SIZE];
  //int bondsMatrix[200][200];
  int dummyInt;
  char fileName[200];
  //  double shift[3]={0.0}, currentCoord;

  dummyInt = sprintf(fileName,"%scoordinates%d_%d.pdb",VMDdir,k,version);
  coordsPtr = fopen(fileName,"w");
  for(i=0;i<N;i++){
    switch(particleState[i]){
    case 4:
      //    if(particleState[i]==4)
      fprintf(coordsPtr,"ATOM      1  C   ALA     1    %8.2f%8.3f%8.2f0.00  1.00      MAIN\n",x[0][i]*scaleFactor,x[1][i]*scaleFactor,x[2][i]*scaleFactor);
      break;
      //    else{
    case 3:
      //if(particleState[i]==3)
      fprintf(coordsPtr,"ATOM      1  C   ARG     1    %8.2f%8.2f%8.2f0.00  1.00      MAIN\n",x[0][i]*scaleFactor,x[1][i]*scaleFactor,x[2][i]*scaleFactor);
      break;
    case 2:
      fprintf(coordsPtr,"ATOM      1  C   ASN     1    %8.2f%8.2f%8.2f0.00  1.00      MAIN\n",x[0][i]*scaleFactor,x[1][i]*scaleFactor,x[2][i]*scaleFactor);
      break;
    case 1:
      fprintf(coordsPtr,"ATOM      1  C   ASP     1    %8.2f%8.2f%8.2f0.00  1.00      MAIN\n",x[0][i]*scaleFactor,x[1][i]*scaleFactor,x[2][i]*scaleFactor);
      break;
    default:
      //     else{	
      fprintf(coordsPtr,"ATOM      1  C   CYS     1    %8.2f%8.2f%8.2f0.00  1.00      MAIN\n",x[0][i]*scaleFactor,x[1][i]*scaleFactor,x[2][i]*scaleFactor);
      break;
    }
  }
  fprintf(coordsPtr,"END");
  fclose(coordsPtr);

  for(i=0;i<N;i++){
    for(j=0;j<N;j++)
      bondsMatrix[i][j] = 0;
  }
  for(i=0;i<N-1;i++)
    bondsMatrix[i][i+1] = 1;

  dummyInt = sprintf(fileName,"%sbonds%d_%d.psf",VMDdir,k,version);
  bondsPtr = fopen(fileName,"w");

  fprintf(bondsPtr,"PSF\n");
  fprintf(bondsPtr,"\n");
  fprintf(bondsPtr,"      11 !NTITLE\n");
  fprintf(bondsPtr,"\n");
  fprintf(bondsPtr,"      %d !NATOM\n",N);
  for(i=1;i<=N;i++){
    switch(particleState[i-1]){
      case 4:
	if(i<10)
	  fprintf(bondsPtr,"       %d MAIN 1    ALA  C    C   0.000000       15.0350           0\n",i);
	else
	  fprintf(bondsPtr,"      %d MAIN 1    ALA  C    C      0.000000       15.0350           0\n",i);
	break;
      case 3:
	if(i<10)
	  fprintf(bondsPtr,"       %d MAIN 1    ARG  C    C   0.000000       15.0350           0\n",i);
	else
	  fprintf(bondsPtr,"      %d MAIN 1    ARG  C    C      0.000000       15.0350           0\n",i);
	break;
      case 2:
	if(i<10)
	  fprintf(bondsPtr,"       %d MAIN 1    ASN  C    C   0.000000       15.0350           0\n",i);
	else
	  fprintf(bondsPtr,"      %d MAIN 1    ASN  C    C      0.000000       15.0350           0\n",i);
	break;
      case 1:
	if(i<10)
	  fprintf(bondsPtr,"       %d MAIN 1    ASP  C    C   0.000000       15.0350           0\n",i);
	else
	  fprintf(bondsPtr,"      %d MAIN 1    ASP  C    C      0.000000       15.0350           0\n",i);
	break;
      default:
	if(i<10)
	  fprintf(bondsPtr,"       %d MAIN 1    CYS  C    C   0.000000       15.0350           0\n",i);
	else
	  fprintf(bondsPtr,"      %d MAIN 1    CYS  C    C      0.000000       15.0350           0\n",i);
	break;
    }
  }
  fprintf(bondsPtr,"\n\n");
  if(totalBonds<10)
    fprintf(bondsPtr,"       %d !NBOND: bonds\n", totalBonds);
  else{
    if(totalBonds<100)
      fprintf(bondsPtr,"      %d !NBOND: bonds\n", totalBonds);
    else
      fprintf(bondsPtr,"     %d !NBOND: bonds\n", totalBonds);
  }
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if(bondsMatrix[i][j] == 1){
	if(i+1<10)
	  fprintf(bondsPtr,"       %d",i+1);
	else{
	  if(i+1<100)
	    fprintf(bondsPtr,"      %d",i+1);
	  else
	    fprintf(bondsPtr,"     %d",i+1);
	}
	if(j+1<10)
	  fprintf(bondsPtr,"       %d",j+1);
	else{
	  if(j+1<100)
	    fprintf(bondsPtr,"      %d",j+1);
	  else
	    fprintf(bondsPtr,"     %d",j+1);
	}
	bondCount++;
	if(bondCount == 4){
	  fprintf(bondsPtr,"\n");
	  bondCount = 0;
	}
      }
    }
  }
  fprintf(bondsPtr,"\n");
  fclose(bondsPtr);

  sprintf(fileName,"%sState%d_%d.vmd",VMDdir,k,version);
  vmdPtr = fopen(fileName,"w");
  fprintf(vmdPtr,"mol delete all\n");
  fprintf(vmdPtr,"mol new bonds%d_%d.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n",k,version);
  fprintf(vmdPtr,"mol addfile coordinates%d_%d.pdb type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n",k,version);
  fprintf(vmdPtr,"mol modstyle 0 all CPK\n");
  fprintf(vmdPtr,"mol modcolor 0 all ResName\n");
  fclose(vmdPtr);

}
