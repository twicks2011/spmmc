#include <math.h>

void hcpGen(double x[][500],int N)
{

	int i,j,k;
	double R = 0.5;

  	/* hcp generator */
  	int m=0;
  	for(i=1;i<=5;i++) {
    	for(j=1;j<=5;j++) {
      		for(k=1;k<=5;k++) {
      			if(i%2==1){
      				if(j%2==1)
      					x[0][m]=2.0*R*k;
      				else
      					x[0][m]=2.0*R*(6-k)-R;
      				x[1][m]=R+(j-1)*sqrt(3)*R;
      			}
      			else{
      				if(j%2==1)
      					x[0][m]=2.0*R*(6-k)-R;
      				else
      					x[0][m]=2.0*R*k;
      				x[1][m]=R+R*sqrt(3.0)/3.0+(5-j)*R*sqrt(3.0);      				
      			}
      			x[2][m]=R+(i-1)*R*2*sqrt(6.0)/3.0;
      			m++;      			
      		}
    	}
  	}
}