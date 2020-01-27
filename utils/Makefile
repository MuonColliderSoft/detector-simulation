CC = g++
LIBS = -llcio -lG4global -lG4physicslists -lG4run -lG4particles -lCLHEP -lCore -lMathCore -lPhysics
IDIR := -I${ILCSOFT}/lcio/HEAD/include \
-I${ILCSOFT}/geant4/10.03.p02/include \
-I${ILCSOFT}/CLHEP/2.3.4.3/include \
-I${ROOTSYS}/include
LDIR := -L${ILCSOFT}/lcio/HEAD/lib \
-L${ILCSOFT}/geant4/10.03.p02/lib64 \
-L${ILCSOFT}/CLHEP/2.3.4.3/lib \
-L${ROOTSYS}/lib


all: mars_to_slcio

mars_to_slcio: mars_to_slcio.cpp
	$(CC) -o $@ $^ $(IDIR) $(LDIR) $(LIBS)
