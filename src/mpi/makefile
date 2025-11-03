#--------------------------------------------------------------
#          Makefile for MPI module
#--------------------------------------------------------------

# List module object filenames
mpi_objects = \
 data.o \
 decomposition.o \
 LLGHeun-mpi.o \
 LLGMidpoint-mpi.o \
 llg_quantum_mpi.o \
 mpi_generic.o \
 mpi_comms.o \
 parallel_rng_seed.o \
 wrapper.o \
 lsf_mpi.o \
 lsf_rk4_mpi.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/mpi/,$(mpi_objects))
