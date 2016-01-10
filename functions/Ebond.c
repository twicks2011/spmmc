#include <math.h>

double Ebond(double *bondAngles, int N){
	int i;
	double Ebond = 0.0;

	double k_theta = 535.2665;
	double theta_0 = 1.91; 

	for(i=1;i<N-1;i++)
		Ebond += k_theta*(bondAngles[i]-theta_0)*(bondAngles[i]-theta_0);

	return Ebond;

}