# simple make file (Migrate.c)
# September 2004
#HOME= /Users/rsgraham
NR=$(HOME)/recipes_c-ansi
#GSL=$(HOME)/GNU_Library


VPATH = ./ : $(NR)/misc/ :  $(NR)/recipes/ : functions/ 


SOURCES= main.c nrutil.c plgndr.c VMDprint.c orderParam.c computeE_torsion.c E_swSum.c E_LJsum.c   initialise.c hcpGen.c rotateAboutU.c bondLengthCheck.c E_LJcrystAv.c E_LJnonCrystAv.c crankShaft.c reptation.c SNorderParam.c Ebond.c positions2bonds.c endRotate.c endBridge.c mersenneTwister.c singleSquareWell.c forwardRept_DE.c backwardRept_DE.c singleParticleMove_DE.c pivotMove_DE.c Estretch.c


OBJECTS= main.o nrutil.o plgndr.o VMDprint.o orderParam.o computeE_torsion.o E_swSum.o E_LJsum.o  initialise.o hcpGen.o rotateAboutU.o bondLengthCheck.o E_LJcrystAv.o E_LJnonCrystAv.o  crankShaft.o reptation.o SNorderParam.o Ebond.o positions2bonds.o endRotate.o endBridge.o mersenneTwister.o singleSquareWell.o forwardRept_DE.o backwardRept_DE.o singleParticleMove_DE.o pivotMove_DE.o Estretch.o

PRODUCT=singleChainOMP

#-------------------------------------------------------------------#
#### icpc compilation ####
CC=icpc
CXX=icpc

CPPFLAGS=-I$(NR)/include -O3 -openmp -ansi #-I$(GSL)/include #-L$(GSL)/lib -O2 #-fast 



#### gcc compilation ####
# To compile with gcc, comment lines 20, 21 and 23 above and uncomment lines 30, 31 and 33 below.

#CC=gcc
#CXX=g++

#CPPFLAGS=-I$(NR)/include -O3 -ansi -fopenmp  #-I$(GSL)/include #-L$(GSL)/lib -O2 #-fast 

#-------------------------------------------------------------------#



all: $(PRODUCT)
$(PRODUCT): $(OBJECTS)
	$(CXX) $(CPPFLAGS) -o $(PRODUCT) $(OBJECTS) -lm 

.c.o:	
	$(CXX)  $(CPPFLAGS)-c $< 

clean:
	rm -f *.o
