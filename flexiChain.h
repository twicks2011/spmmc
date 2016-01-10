#ifndef _SINGLE_CHAIN_H_
#define _SINGLE_CHAIN_H_

#include <stdio.h>

void initialise(double[][500], int, char*);
void hcpGen(double[][500],int);
void positions2bonds(double[][500], double*, double*, int);
void rotateAboutU(double[][500], int, int, double);
void crankShaft(double[][500], int, int, double);
void reptation(double[][500], int, int, double, double); 
void endRotate(double[][500], int, int, double, double, double);
double endBridge(double[][500], int, int, double, double, int *);
double E_LJsum(double[][500], int);
double E_swSum(double[][500], int, double);
double Ebond(double*, int);
double E_LJpartial(double **, int, int);
double E_LJsubchain(double**, int, int);
double E_LJcrystAv(double[][500], int, int*);
double E_LJnonCrystAv(double[][500], int, int*);
double computeE_torsion(double*, int, int, int);
double Estretch(double[][500], int, double, int, int);
void VMDprint(double[][500], int, int, int, int*, char*);
int orderParam(double[][500], int, int, int*);
int SNorderParam(double[][500], int, double, int*,double);
int bondLengthCheck(double[][500], int);
double E_swSingle(double rijsq,  double lambdaSq);
double forwardRep_DE(double x[][500],double x_old[][500], int N, double lambdaSq);
double backwardRep_DE(double x[][500],double x_old[][500], int N, double lambdaSq);
double singleParticleMove_DE(double x[][500],double x_old[][500], int N, int movedParticle, double lambdaSq);
double pivotMove_DE(double x[][500],double x_old[][500], int N, int pivotParticle, double lambdaSq);
//int DETorderParam(double[][500], int, int, int, double, int, int, int*);

#endif
