#include <math.h>
#include <stdio.h>
#include "../flexiChain.h"

double E_swSum(double x[][500], int N, double lambda)
{
	double Ei, Esum=0.0;
	int i, j, cmpt;
	double xi[3], xij[3], rijsq;

	double lambdaSQ = lambda*lambda;

	for(i=0;i<N-1;i++){
		Ei = 0.0;
		for(cmpt=0;cmpt<3;cmpt++)
			xi[cmpt] = x[cmpt][i];
		for(j=i+1;j<N;j++){
			rijsq = 0.0;
			for(cmpt=0;cmpt<3;cmpt++){
				xij[cmpt] = x[cmpt][j]-xi[cmpt];
				rijsq += xij[cmpt]*xij[cmpt];
			}
			//printf("rijsq[%d][%d]: %f\n", i, j, rijsq);			
			if((i-j)*(i-j)>1){

			   Ei    +=   E_swSingle( rijsq  ,  lambdaSQ) ;
			
			}
		}
		//printf("Ei[%d]: %f\n", i, Ei);
		Esum += Ei;
	}

	return Esum;
}
