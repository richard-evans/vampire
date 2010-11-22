#===================================================================
#
#								Makefile for VAMPIRE
#
#===================================================================
#export OMPI_CXX=icc
#export OMPI_CXX=pathCC
# Compilers
ICC=icc -DCOMP='"Intel C++ Compiler"'
GCC=g++ -DCOMP='"GNU C++ Compiler"'
PCC=pathCC -DCOMP='"Pathscale C++ Compiler"'
MPICC=mpic++ -DMPICF
export LANG=C
export LC_ALL=C

# Debug Flags
ICC_DBCFLAGS= -O0 -C -I./hdr
ICC_DBLFLAGS= -lstdc++ -C -I./hdr

GCC_DBCFLAGS= -Wall -Wextra -O0 -fbounds-check -pedantic -ansi -Wno-long-long -I./hdr
GCC_DBLFLAGS= -lstdc++ -fbounds-check -I./hdr

PCC_DBCFLAGS= -O0 -I./hdr
PCC_DBLFLAGS= -O0 -I./hdr

# Performance Flags
ICC_CFLAGS= -O3 -axSSE3 -ipo -static -fno-alias -align -falign-functions -I./hdr
ICC_LDFLAGS= -lstdc++ -ipo -I./hdr -axSSE3
#ICC_CFLAGS= -O3 -xT -ipo -static -fno-alias -align -falign-functions -vec-report -I./hdr
#ICC_LDFLAGS= -lstdc++ -ipo -I./hdr -xT -vec-report

GCC_CFLAGS=-O3 -msse3 -falign-labels -falign-loops -funroll-all-loops -fexpensive-optimizations -funroll-loops -I./hdr
GCC_LDFLAGS= -lstdc++ -I./hdr

PCC_CFLAGS=-O2 -march=barcelona -ipa -I./hdr
PCC_LDFLAGS= -lstdc++ -I./hdr -O2 -march=barcelona -ipa

# Objects
OBJECTS= \
obj/create/create_system2.o \
obj/create/cs_create_crystal_structure2.o \
obj/create/cs_create_system_type2.o \
obj/create/cs_create_neighbour_list2.o \
obj/create/cs_particle_shapes.o \
obj/create/cs_set_atom_vars2.o \
obj/create/cs_voronoi2.o \
obj/data/atoms.o \
obj/data/grains.o \
obj/main/initialise_variables.o \
obj/main/main.o \
obj/main/material.o \
obj/mpi/LLG_mpi.o \
obj/mpi/mpi_generic.o \
obj/mpi/mpi_create2.o \
obj/mpi/mpi_comms.o \
obj/program/bmark.o \
obj/program/curie_temperature.o \
obj/program/hysteresis.o \
obj/program/LLB_Boltzmann.o \
obj/program/static_hysteresis.o \
obj/random/mtrand.o \
obj/random/random.o \
obj/simulate/fields.o \
obj/simulate/demag.o \
obj/simulate/LLB.o \
obj/simulate/LLG.o \
obj/simulate/simulate_system.o \
obj/simulate/standard_programs.o \
obj/utility/errors.o \
obj/utility/statistics.o \
obj/utility/vio.o \
obj/utility/vout.o \
obj/utility/vmath.o \
obj/utility/units.o 

ICC_OBJECTS=$(OBJECTS:.o=_i.o)
ICCDB_OBJECTS=$(OBJECTS:.o=_idb.o)
GCCDB_OBJECTS=$(OBJECTS:.o=_gdb.o)
PCCDB_OBJECTS=$(OBJECTS:.o=_pdb.o)

MPI_OBJECTS=$(OBJECTS:.o=_mpi.o)
MPI_ICC_OBJECTS=$(OBJECTS:.o=_i_mpi.o)
MPI_PCC_OBJECTS=$(OBJECTS:.o=_p_mpi.o)
MPI_ICCDB_OBJECTS=$(OBJECTS:.o=_idb_mpi.o)
MPI_GCCDB_OBJECTS=$(OBJECTS:.o=_gdb_mpi.o)
MPI_PCCDB_OBJECTS=$(OBJECTS:.o=_pdb_mpi.o)

EXECUTABLE=vampire

all: $(OBJECTS) gcc

# Serial Targets
gcc: $(OBJECTS)
	$(GCC) $(GCC_LDFLAGS) $(OBJECTS) -o $(EXECUTABLE)

$(OBJECTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $<

intel: $(ICC_OBJECTS)
	$(ICC) $(ICC_LDFLAGS) $(ICC_OBJECTS) -o $(EXECUTABLE)

$(ICC_OBJECTS): obj/%_i.o: src/%.cpp
	$(ICC) -c -o $@ $(ICC_CFLAGS) $<

gcc-debug: $(GCCDB_OBJECTS)
	$(GCC) $(GCC_DBLFLAGS) $(GCCDB_OBJECTS) -o $(EXECUTABLE)

$(GCCDB_OBJECTS): obj/%_gdb.o: src/%.cpp
	$(GCC) -c -o $@ $(GCC_DBCFLAGS) $<

intel-debug: $(ICCDB_OBJECTS)
	$(ICC) $(ICC_DBLFLAGS) $(ICCDB_OBJECTS) -o $(EXECUTABLE)

$(ICCDB_OBJECTS): obj/%_idb.o: src/%.cpp
	$(ICC) -c -o $@ $(ICC_DBCFLAGS) $<

pathscale-debug: $(ICCDB_OBJECTS)
	$(PCC) $(PCC_DBLFLAGS) $(PCCDB_OBJECTS) -o $(EXECUTABLE)

$(PCCDB_OBJECTS): obj/%_pdb.o: src/%.cpp
	$(PCC) -c -o $@ $(PCC_DBCFLAGS) $<

# MPI Targets

mpi-gcc: $(MPI_OBJECTS)
	$(MPICC) $(GCC_LDFLAGS) $(MPI_OBJECTS) -o $(EXECUTABLE)
#export OMPI_CXX=icc
$(MPI_OBJECTS): obj/%_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(GCC_CFLAGS) $<

mpi-intel: $(MPI_ICC_OBJECTS)
	$(MPICC) $(ICC_LDFLAGS) $(MPI_ICC_OBJECTS) -o $(EXECUTABLE)
$(MPI_ICC_OBJECTS): obj/%_i_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(ICC_CFLAGS) $<

mpi-pathscale: $(MPI_PCC_OBJECTS)
	$(MPICC) $(PCC_LDFLAGS) $(MPI_PCC_OBJECTS) -o $(EXECUTABLE)
$(MPI_PCC_OBJECTS): obj/%_p_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(PCC_CFLAGS) $<

mpi-gcc-debug: $(MPI_GCCDB_OBJECTS)
	$(MPICC) $(GCC_DBLFLAGS) $(MPI_GCCDB_OBJECTS) -o $(EXECUTABLE)

$(MPI_GCCDB_OBJECTS): obj/%_gdb_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(GCC_DBCFLAGS) $<

mpi-intel-debug: $(MPI_ICCDB_OBJECTS)
	$(MPICC) $(ICC_DBLFLAGS) $(MPI_ICCDB_OBJECTS) -o $(EXECUTABLE)

$(MPI_ICCDB_OBJECTS): obj/%_idb_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(ICC_DBCFLAGS) $<

mpi-pathscale-debug: $(MPI_PCCDB_OBJECTS)
	$(MPICC) $(PCC_DBLFLAGS) $(MPI_PCCDB_OBJECTS) -o $(EXECUTABLE)

$(MPI_PCCDB_OBJECTS): obj/%_pdb_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(PCC_DBCFLAGS) $<

clean:
	@rm -f obj/*.o
	@rm -f obj/*/*.o

purge:
	@rm -f obj/*.o
	@rm -f obj/*/*.o
	@rm -f vampire

tidy:	
	@rm -f *~
	@rm -f src/*~
	@rm -f src/*/*~
	@rm -f hdr/*~

package:
	@bash .pack


