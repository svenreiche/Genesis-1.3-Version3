#
# GENESIS makefile
#
# libraries
#
#LIB= -lmpi_cxx -lgfortran -lmpe -L/gpfs/home/reiche/Apps/MPE/lib -llmpe -lmpe
LIB= -lm -lmpi_cxx -lz -lgfortran
#
#INCLUDE=-I./include -I/gpfs/home/reiche/Apps/MPE/include
INCLUDE=-I./include
#
# cpp-macros
#
DMACRO = -DHARMMAX=7 -DNZMAX=10000 -DGENESIS_VERSION=3.0 -DMPEEEE
#
#  compilers
#
VPATH = src/C++ src/C++/Input src/C++/Output src/Fortran src/C++/Util
CCOMPILER = h5pcc
FCOMPILER = gfortran
#
#  executable name
#
EXECUTABLE = genesis
#
# targets
#
OBJECTS = check.o diagno.o esource.o field.o incoherent.o \
 math.o partsim.o pushp.o loadbeam.o loadrad.o magfield.o \
 tdepend.o track.o string.o rpos.o scan.o source.o stepz.o importbeam.o\
 timerec.o initrun.o  input.o output.o mpi_common.o fort2cpp.o \
 main.o outputbase.o  outputHDF5.o outdumpbase.o  \
 outdumpHDF5.o readparticlebase.o readparticleHDF5.o HDF5base.o \
 readfieldbase.o readfieldHDF5.o RandomU.o


genesis:	$(OBJECTS)
	$(CCOMPILER)  -o $(EXECUTABLE) $(OBJECTS) $(LIB) $(MPELIB)



.cpp.o:
	$(CCOMPILER) -O -c $(DMACRO) $(INCLUDE) $(MPEINCLUDE) $<

.f.o:
	$(FCOMPILER) -O  -c $<

clean:
	rm -f *~
	rm -f *.o
	rm $(EXECUTABLE)

install:
	cp ./$(EXECUTABLE) ../bin/

beta:
	cp ./$(EXECUTABLE) ~/bin/$(EXECUTABLE).beta











