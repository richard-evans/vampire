#--------------------------------------------------------------
#          Makefile for main module
#--------------------------------------------------------------

# List module object filenames
main_objects =\
command.o \
githash.o \
initialise_variables.o \
main.o \
material.o \
version.o

# Append module objects to global tree
OBJECTS+=$(addprefix obj/main/,$(main_objects))

# Add specific options for every compiler (ugh)
# GNU
obj/main/githash.o : OPTIONS = -DGHASH='$(GHASH)'
obj/main/githash_mpi.o : OPTIONS = -DGHASH='$(GHASH)'
# GNU DEBUG
obj/main/githash_gdb.o : OPTIONS = -DGHASH='$(GHASH)'
obj/main/githash_gdb_mpi.o : OPTIONS = -DGHASH='$(GHASH)'
# LLVM
obj/main/githash_llvm.o : OPTIONS = -DGHASH='$(GHASH)'
obj/main/githash_llvm_mpi.o : OPTIONS = -DGHASH='$(GHASH)'
# INTEL
obj/main/githash_i.o : OPTIONS = -DGHASH='$(GHASH)'
obj/main/githash_i_mpi.o : OPTIONS = -DGHASH='$(GHASH)'
# ARCHER
obj/main/githash_archer_mpi.o : OPTIONS = -DGHASH='$(GHASH)'
