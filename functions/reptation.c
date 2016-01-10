#include <math.h>
#include <time.h>

#define PI 3.141592653589793238462643383279502884197

void reptation(double x[][500], int N, int forward, double rndnum_1, double rndnum_2)
{	
	int i, j, cmpt;
	double r[3];	
	double rndnum_1sq = rndnum_1*rndnum_1;
	double rndnum_2sq = rndnum_2*rndnum_2;
	
	r[0] = 2.0*rndnum_1*sqrt(1.0-rndnum_1sq-rndnum_2sq);
	r[1] = 2.0*rndnum_2*sqrt(1.0-rndnum_1sq-rndnum_2sq);
	r[2] = 1.0-2.0*(rndnum_1sq+rndnum_2sq);

	if(forward == 1){
		/*FORWARD REPTATION*/
		for(i=1;i<N;i++){
			for(cmpt=0;cmpt<3;cmpt++)
				x[cmpt][i-1] = x[cmpt][i];
		}
		for(cmpt=0;cmpt<3;cmpt++)
			x[cmpt][N-1] = x[cmpt][N-2] + r[cmpt];
	}
	else{
		/*BACKWARD REPTATION*/
		for(i=N-1;i>0;i--){
			for(cmpt=0;cmpt<3;cmpt++)
				x[cmpt][i] = x[cmpt][i-1];

		}
		for(cmpt=0;cmpt<3;cmpt++)
			x[cmpt][0] = x[cmpt][1] + r[cmpt];
	}

}
