#--------------------------------------------------------------
#          Makefile for dipole module
#--------------------------------------------------------------

# List module object filenames
dipole_objects =\
atomistic.o \
data.o \
energy.o \
field.o \
info.o \
get_cells_properties.o \
get_tensor.o \
initialize.o \
interface.o \
inter.o \
intra.o \
macrocell.o \
memory.o \
mpi.o \
mpi2.o \
output_atomistic_field.o \
tensor.o \
update.o \
fft_macrocell.o \
fft_atomistic.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/dipole/,$(dipole_objects))
