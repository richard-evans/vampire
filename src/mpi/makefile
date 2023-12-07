#--------------------------------------------------------------
#          Makefile for MPI module
#--------------------------------------------------------------

# List module object filenames
mpi_objects =\
data.o \
decomposition.o \
LLGHeun-mpi.o \
LLGMidpoint-mpi.o \
mpi_generic.o \
mpi_comms.o \
parallel_rng_seed.o \
wrapper.o \
LSF-mpi.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/mpi/,$(mpi_objects))
