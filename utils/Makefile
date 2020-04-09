CC = g++
LIBS = -llcio -lG4global -lG4physicslists -lG4run -lG4particles -lCLHEP -lCore -lMathCore -lPhysics
IDIR := -I${LCIO}/include \
-I${CLHEP_INCLUDE_DIR} \
-I${G4INSTALL}/include \
-I$(shell root-config --incdir)
LDIR := -L${LCIO}/lib \
-L${CLHEP_DIR}/lib \
-L${G4INSTALL}/lib64 \
-L$(shell root-config --libdir)


all: mars_to_slcio

mars_to_slcio: mars_to_slcio.cpp
	$(CC) -o $@ $^ $(IDIR) $(LDIR) $(LIBS)
