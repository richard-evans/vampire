#--------------------------------------------------------------
#          Makefile for montecarlo module
#--------------------------------------------------------------

# List module object filenames
montecarlo_objects =\
data.o \
initialize.o \
interface.o \
mc.o \
mc_moves.o \
cmc.o \
masked_cmc_mc.o \
cmc_mc.o \
monte_carlo_preconditioning.o \
mc-mpi.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/montecarlo/,$(montecarlo_objects))
