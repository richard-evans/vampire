#--------------------------------------------------------------
#          Makefile for voro++ module
#--------------------------------------------------------------

# List module object filenames
voro++_objects =\
cell.o \
common.o \
container.o \
unitcell.o \
v_compute.o \
c_loops.o \
v_base.o \
wall.o \
pre_container.o \
container_prd.o

include src/voro++/Makefile.dep

# Append module objects to global tree
OBJECTS+=$(addprefix obj/voro++/,$(voro++_objects))
