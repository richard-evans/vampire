#--------------------------------------------------------------
#          Makefile for dipole module
#--------------------------------------------------------------

# List module object filenames
dipole_objects =\
atomistic.o \
data.o \
energy.o \
field.o \
fft_update.o \
initialize.o \
interface.o \
inter.o \
intra.o \
macrocell.o \
memory.o \
mpi.o \
output_atomistic_field.o \
tensor.o \
update.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/dipole/,$(dipole_objects))
