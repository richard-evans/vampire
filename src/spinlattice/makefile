#--------------------------------------------------------------
#          Makefile for sld module
#--------------------------------------------------------------

# List module object filenames
spinlattice_objects =\
data.o \
initialize.o \
interface.o\
forces.o\
fields.o\
suzuki_trotter.o\
temperatures.o\
energy.o\
suzuki_trotter_mpi.o\
tests.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/spinlattice/,$(spinlattice_objects))
