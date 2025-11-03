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
suzuki-trotter.o\
temperatures.o\
energy.o\
suzuki-trotter-mpi.o\
tests.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/spinlattice/,$(spinlattice_objects))
