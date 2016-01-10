#include <math.h>
#include <stdio.h>
//#include <nrutil.h>
//#include <nr_RSG.h>

#define L 6
#define DC 0.5

double plgndr(int, int, double);

int SNorderParam(double x[][500], int N, double lambda, int *particleState, double E_sw){

	if(E_sw < 0.0)
		return -(int)E_sw;
	else
		return 0;

	int crystParts = 0;
	int i,j,k,m,cmpt;
	double xij[3], rijsq;
	int Nn[200]={0}, Np[200][200]={0};
	double Sct[200][200]={0.0}, Sp[200][200]={0.0};
	double qr[200][200]={0.0}, qi[200][200]={0.0};
	double d[200][200]={0.0};
	int Nc[200]={0};
	double mqi2, mqj2,qiqjc;
	int cryst[200]={0}, coil[200]={0};

	double lambdaSQ = lambda*lambda;	
	double c;
	double factRatio = 1.0, sign;

/*
	return floor(  sqrt((x[0][0]-x[0][N-1])*(x[0][0]-x[0][N-1]) +
			    (x[1][0]-x[1][N-1])*(x[1][0]-x[1][N-1]) +
			    (x[2][0]-x[2][N-1])*(x[2][0]-x[2][N-1]) ));
*/

	
	for(i=0;i<N;i++){
		particleState[i] = 2;
		k=0;
		for(j=0;j<N;j++){
			if(j != i){
				rijsq = 0.0;
				for(cmpt=0;cmpt<3;cmpt++){
					xij[cmpt] = x[cmpt][j]-x[cmpt][i];
					rijsq += xij[cmpt]*xij[cmpt];
				}
				if(rijsq < lambdaSQ){					
					Nn[i]++;
					Np[i][k]=j;
					Sct[i][k]=xij[2]/(sqrt(rijsq));
					Sp[i][k]=atan2(xij[1],xij[0]);
					k++;
				}
			}
		}
	}
	/*
	for(i=0;i<N;i++){	
		if(Nn[i] >= 5){
			particleState[i] = 1;		
			crystParts++;
		}
	}	
	return crystParts;
	*/

	double Lfactor = sqrt((2.0*L+1.0)/(4*M_PI));
	for(i=0;i<N;i++){
		sign = 1.0;
		for(m=0;m<=L;m++){
			if(m>0){
				factRatio *= (L+m)*(L-m+1);
				sign *= -1.0;
			}			
			for(j=0;j<Nn[i];j++){
				c = Lfactor/sqrt(factRatio);
				c *= plgndr(L,m,Sct[i][j]);
				qr[i][m+L] += c*sin(m*Sp[i][j])/(1.0*Nn[i]);
				qi[i][m+L] += c*sin(m*Sp[i][j])/(1.0*Nn[i]);
				if(m>0){
					qr[i][L-m] += sign*c*cos(m*Sp[i][j])/(1.0*Nn[i]);
					qi[i][L-m] -= sign*c*cos(m*Sp[i][j])/(1.0*Nn[i]);
				}
			}
		}
		factRatio = 1.0;
	}

	for(i=0;i<N;i++){
		for(j=0;j<Nn[i];j++){
			mqi2 = 0.0;
			mqj2 = 0.0;
			qiqjc = 0.0;
			for(m=-L;m<=L;m++){
			        mqi2 += qr[i][m+L]*qr[i][m+L]+qi[i][m+L]*qi[i][m+L];
				mqj2 += qr[Np[i][j]][m+L]*qr[Np[i][j]][m+L]+qi[Np[i][j]][m+L]*qi[Np[i][j]][m+L];
 				qiqjc += (qr[i][m+L]*qr[Np[i][j]][m+L]) + (qi[i][m+L]*qi[Np[i][j]][m+L]);
			}
			d[i][j]=qiqjc/sqrt(mqi2*mqj2);
		}
	}

	for(i=0;i<N;i++){
		for(j=0;j<Nn[i];j++){
			if(d[i][j]-DC>0)
				Nc[i]++;
		}
		if(Nn[i] >= 5 && Nc[i] >= Nn[i]-1){			
			cryst[i] = 1;
			particleState[i] = 1; //crystalline particle
			crystParts++;
		}
		else{
			if(Nn[i] <= 4){
				coil[i] = 1; 
				particleState[i] = 0; //coil-like
			}
			else
				particleState[i] = 4; //intermediate						
		}
	}

	return crystParts;
}

double plgndr(int l, int m, double x)
{
  //void nrerror(char error_text[]);
  double fact,pll,pmm,pmmp1,somx2;
  int i,ll;

  //if (m < 0 || m > l || fabs(x) > 1.0)
  //  nrerror("Bad arguments in routine plgndr");
  pmm=1.0;
  if (m > 0) {
    somx2=sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for (i=1;i<=m;i++) {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else {
    pmmp1=x*(2*m+1)*pmm;
    if (l == (m+1))
      return pmmp1;
    else {
      for (ll=m+2;ll<=l;ll++) {
        pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
        pmm=pmmp1;
        pmmp1=pll;
      }
      return pll;
    }
  }
}
