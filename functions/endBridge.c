#include <stdio.h>
#include <math.h>

double endBridge(double x[][500], int N, int end, double rndnum_1, double rndnum_2, int *changed){
	int i, j, cmpt, END;

	double x_old[3][500];

	//store original configuration
	for(cmpt=0;cmpt<3;cmpt++){
		for(i=0;i<N;i++)
			x_old[cmpt][i] = x[cmpt][i];
	}


	if(end==0)
		END = 0;
	else
		END=N-1;

	double rijsq;
	int numNeighbours=0, neighbourIndices[500];


	//locate and count potential bridging sites
	for(j=3;j<N-3;j++){		
		rijsq=0.0;
		for(cmpt=0;cmpt<3;cmpt++)
			rijsq += (x[cmpt][END]-x[cmpt][j])*(x[cmpt][END]-x[cmpt][j]);
		if(rijsq < 4.0){			
			neighbourIndices[numNeighbours]=j;
			numNeighbours++;
		}		
	}
	
	int b_a = numNeighbours;

	if(b_a==0){
	  *changed = N;
		return 0.0;	
	}


	//========Choosing and relabelling====
	//choose the particle to remove
	int removedPart = neighbourIndices[(int)(rndnum_1*numNeighbours)];
	double R_a=0.0, lTemp;

	//compute distance from end to chosen particle
	for(cmpt=0;cmpt<3;cmpt++){
		lTemp = x[cmpt][removedPart]-x[cmpt][END];
		R_a += lTemp*lTemp;			
	}
	R_a = sqrt(R_a);

	
	if(end==0){				
		removedPart--;
		//relabel all particles up to removed one
		for(j=0;j<removedPart;j++){
			for(cmpt=0;cmpt<3;cmpt++)
				x[cmpt][j] = x_old[cmpt][removedPart-1-j];		
		}
	}
	
	else{				
		removedPart++;
		for(j=removedPart+1;j<N;j++){
			for(cmpt=0;cmpt<3;cmpt++)
				x[cmpt][j] = x_old[cmpt][N+removedPart-j];			
		}
	}	


	
	//========Reinsert the chosen particle====
	double l[3];
	double midPoint[3];
	//find midpoint between 
	for(cmpt=0;cmpt<3;cmpt++){
		l[cmpt] = x[cmpt][removedPart+1]-x[cmpt][removedPart-1];
		midPoint[cmpt] = 0.5*(x[cmpt][removedPart-1]+x[cmpt][removedPart+1]);		
	}
	double modLsq = l[0]*l[0]+l[1]*l[1]+l[2]*l[2];
	double r[3];

	double lsq = l[1]*l[1] + l[2]*l[2];

	if( lsq > 1e-12){
	  r[0] = 0.0;
	  r[1] = -l[2]/sqrt(lsq);
	  r[2] = l[1]/sqrt(lsq);
	}else{
	  r[0]=0.0;
	  r[1]=1.0;
	  r[2]=0.0;
	}

	//printf("%f %f %f\n",l[0],l[1],l[2]);
	
	for(cmpt=0;cmpt<3;cmpt++){
		r[cmpt] *= sqrt(1-0.25*modLsq);
		x[cmpt][removedPart] = midPoint[cmpt]+r[cmpt];
	}

	
	double phi = (2.0*rndnum_2-1.0)*M_PI;
	double cp = cos(phi);
  	double omcp = 1.0-cp;
  	double sp = sin(phi);
  	double rotationMatrix[3][3];
  	double U[3];
  	for(cmpt=0;cmpt<3;cmpt++)
  		U[cmpt] = l[cmpt]/sqrt(modLsq);

  	rotationMatrix[0][0] = cp + U[0]*U[0]*omcp;
  	rotationMatrix[0][1] = U[0]*U[1]*omcp-U[2]*sp;
  	rotationMatrix[0][2] = U[0]*U[2]*omcp+U[1]*sp;
  	rotationMatrix[1][0] = U[0]*U[1]*omcp+U[2]*sp;
  	rotationMatrix[1][1] = cp+U[1]*U[1]*omcp;
  	rotationMatrix[1][2] = U[1]*U[2]*omcp-U[0]*sp;
  	rotationMatrix[2][0] = U[0]*U[2]*omcp-U[1]*sp;
  	rotationMatrix[2][1] = U[1]*U[2]*omcp+U[0]*sp;
  	rotationMatrix[2][2] = cp+U[2]*U[2]*omcp;

  	double r_bond[3];
  	for(cmpt=0;cmpt<3;cmpt++)
  		r_bond[cmpt] = x[cmpt][removedPart]-x[cmpt][removedPart-1];
  	double x_rotated[3] = {0.0};
  	for(i=0;i<3;i++){
  		for(j=0;j<3;j++)
  			x_rotated[i] += rotationMatrix[i][j]*r_bond[j];
  		x_rotated[i] += x[i][removedPart-1];
  		x[i][removedPart] = x_rotated[i];
  	}
	
  	int b_b = 0;
  	for(j=3;j<N-3;j++){		
		rijsq=0.0;
		for(cmpt=0;cmpt<3;cmpt++)
			rijsq += (x[cmpt][END]-x[cmpt][j])*(x[cmpt][END]-x[cmpt][j]);
		if(rijsq < 4.0){			
			b_b++;			
		}		
	}
	double R_b=0.0;
	if(end==0){
		for(cmpt=0;cmpt<3;cmpt++){
			lTemp = x[cmpt][removedPart+1]-x[cmpt][0];
			R_b += lTemp*lTemp;
		}
		R_b = sqrt(R_b);
	}
	else{
		for(cmpt=0;cmpt<3;cmpt++){
			lTemp = x[cmpt][removedPart-1]-x[cmpt][N-1];
			R_b += lTemp*lTemp;
		}
		R_b = sqrt(R_b);
	}	

	*changed = removedPart;
	
	return log((1.0*b_a/(1.0*b_b))*(R_b/R_a));	
}
