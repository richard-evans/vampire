#--------------------------------------------------------------
#          Makefile for quantum module
#--------------------------------------------------------------

# List module object filenames
quantum_objects =\
data.o \
initialize.o \
interface.o \
noise.o \
llg.o \
llg-FFT-serial.o \
llg-FFT-mpi.o \
llg-HO-serial.o \
llg-HO-mpi.o \
RK4.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/quantum/,$(quantum_objects))
