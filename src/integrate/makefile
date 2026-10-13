#--------------------------------------------------------------
#          Makefile for integrate module
#--------------------------------------------------------------

# List module object filenames
integrate_objects =\
data.o \
increment_time.o \
initialize.o \
interface.o \
integrate.o \
llb.o \
llg_heun.o \
llg_heun_mpi.o \
llg_midpoint.o \
llg_midpoint_mpi.o \
llg_quantum.o \
lsf.o \
lsf_mpi.o \
lsf_rk4.o \
lsf_rk4_mpi.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/integrate/,$(integrate_objects))
