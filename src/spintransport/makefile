#--------------------------------------------------------------
#          Makefile for spintransport module
#--------------------------------------------------------------

# List module object filenames
spintransport_objects =\
data.o \
field.o \
initialize.o \
interface.o \
magnetization.o \
resistance.o \
update.o \
get_old_spin_array.o \
output_atomistic_spin_pumping.o \
spin_pumping.o \
voltage.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/spintransport/,$(spintransport_objects))
