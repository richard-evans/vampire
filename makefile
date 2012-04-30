#===================================================================
#
#								Makefile for VAMPIRE
#
#===================================================================

voronoi_path='"qvoronoi"'

export OMPI_CXX=g++
#export OMPI_CXX=icc
#export OMPI_CXX=pathCC
#export MPICH_CXX=g++
export MPICH_CXX=bgxlc++
# Compilers
ICC=icc -DCOMP='"Intel C++ Compiler"' -DVORONOI=$(voronoi_path)
GCC=g++ -DCOMP='"GNU C++ Compiler"' -DVORONOI=$(voronoi_path)
PCC=pathCC -DCOMP='"Pathscale C++ Compiler"' -DVORONOI=$(voronoi_path)
IBM=bgxlc++ -DCOMP='"IBM XLC++ Compiler"' -DVORONOI=$(voronoi_path)
MPICC=mpicxx -DMPICF -DVORONOI=$(voronoi_path)

export LANG=C
export LC_ALL=C

# LIBS
LIBS=-lstdc++
CUDALIBS=-L/usr/local/cuda/lib64/ -lcuda -lcudart
# Debug Flags
ICC_DBCFLAGS= -O0 -C -I./hdr
ICC_DBLFLAGS= -C -I./hdr

GCC_DBCFLAGS= -Wall -Wextra -O0 -fbounds-check -pedantic -std=c++98 -Wno-long-long -I./hdr
GCC_DBLFLAGS= -lstdc++ -fbounds-check -I./hdr

PCC_DBCFLAGS= -O0 -I./hdr
PCC_DBLFLAGS= -O0 -I./hdr

IBM_DBCFLAGS= -O0 -Wall -pedantic -Wextra -I./hdr
IBM_DBLFLAGS= -O0 -Wall -pedantic -Wextra -I./hdr

# Performance Flags
ICC_CFLAGS= -O3 -axSSE3 -fno-alias -align -falign-functions -I./hdr
ICC_LDFLAGS= -I./hdr -axSSE3
#ICC_CFLAGS= -O3 -xT -ipo -static -fno-alias -align -falign-functions -vec-report -I./hdr
#ICC_LDFLAGS= -lstdc++ -ipo -I./hdr -xT -vec-report

GCC_CFLAGS=-O3 -falign-labels -falign-loops -funroll-all-loops -fexpensive-optimizations -funroll-loops -I./hdr
GCC_LDFLAGS= -lstdc++ -I./hdr

PCC_CFLAGS=-O2 -march=barcelona -ipa -I./hdr
PCC_LDFLAGS= -I./hdr -O2 -march=barcelona -ipa

NVCC_FLAGS=-I/usr/local/cuda/include -I./hdr --compiler-bindir=/usr/bin/g++-4.2 --compiler-options=-O3,-DCUDA  --ptxas-options=-v --maxrregcount=32 -arch=sm_13 -O3 
NVCC=nvcc -DCOMP='"GNU C++ Compiler"'

IBM_CFLAGS=-O5 -qarch=450 -qtune=450 -I./hdr 
IBM_LDFLAGS= -lstdc++ -I./hdr -O5 -qarch=450 -qtune=450


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
obj/data/category.o \
obj/data/cells.o \
obj/data/grains.o \
obj/main/initialise_variables.o \
obj/main/main.o \
obj/main/material.o \
obj/mpi/LLGHeun-mpi.o \
obj/mpi/LLGMidpoint-mpi.o \
obj/mpi/mpi_generic.o \
obj/mpi/mpi_create2.o \
obj/mpi/mpi_comms.o \
obj/program/bmark.o \
obj/program/cmc_anisotropy.o \
obj/program/curie_temperature.o \
obj/program/diagnostics.o \
obj/program/field_cool.o \
obj/program/hamr.o \
obj/program/hybrid_cmc.o \
obj/program/hysteresis.o \
obj/program/LLB_Boltzmann.o \
obj/program/static_hysteresis.o \
obj/program/time_series.o \
obj/program/two_temperature.o \
obj/random/mtrand.o \
obj/random/random.o \
obj/simulate/energy.o \
obj/simulate/fields.o \
obj/simulate/demag.o \
obj/simulate/LLB.o \
obj/simulate/LLGHeun.o \
obj/simulate/LLGMidpoint.o \
obj/simulate/mc.o \
obj/simulate/cmc.o \
obj/simulate/cmc_mc.o \
obj/simulate/sim.o \
obj/simulate/standard_programs.o \
obj/utility/errors.o \
obj/utility/statistics.o \
obj/utility/vconfig.o \
obj/utility/vio.o \
obj/utility/vmath.o \
obj/utility/units.o

ICC_OBJECTS=$(OBJECTS:.o=_i.o)
IBM_OBJECTS=$(OBJECTS:.o=_ibm.o)
ICCDB_OBJECTS=$(OBJECTS:.o=_idb.o)
GCCDB_OBJECTS=$(OBJECTS:.o=_gdb.o)
PCCDB_OBJECTS=$(OBJECTS:.o=_pdb.o)
IBMDB_OBJECTS=$(OBJECTS:.o=_ibmdb.o)

MPI_OBJECTS=$(OBJECTS:.o=_mpi.o)
MPI_ICC_OBJECTS=$(OBJECTS:.o=_i_mpi.o)
MPI_PCC_OBJECTS=$(OBJECTS:.o=_p_mpi.o)
MPI_IBM_OBJECTS=$(OBJECTS:.o=_ibm_mpi.o)
MPI_ICCDB_OBJECTS=$(OBJECTS:.o=_idb_mpi.o)
MPI_GCCDB_OBJECTS=$(OBJECTS:.o=_gdb_mpi.o)
MPI_PCCDB_OBJECTS=$(OBJECTS:.o=_pdb_mpi.o)
MPI_IBMDB_OBJECTS=$(OBJECTS:.o=_ibmdb_mpi.o)

CUDA_OBJECTS=$(OBJECTS:.o=_cuda.o)
EXECUTABLE=zspin

all: $(OBJECTS) gcc

# Serial Targets
gcc: $(OBJECTS)
	$(GCC) $(GCC_LDFLAGS) $(LIBS) $(OBJECTS) -o $(EXECUTABLE)

$(OBJECTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $<

intel: $(ICC_OBJECTS)
	$(ICC) $(ICC_LDFLAGS) $(LIBS) $(ICC_OBJECTS) -o $(EXECUTABLE)

$(ICC_OBJECTS): obj/%_i.o: src/%.cpp
	$(ICC) -c -o $@ $(ICC_CFLAGS) $<

ibm: $(IBM_OBJECTS)
	$(IBM) $(IBM_LDFLAGS) $(IBM_OBJECTS) -o $(EXECUTABLE)

$(IBM_OBJECTS): obj/%_ibm.o: src/%.cpp
	$(IBM) -c -o $@ $(IBM_CFLAGS) $<

gcc-debug: $(GCCDB_OBJECTS)
	$(GCC) $(GCC_DBLFLAGS) $(LIBS) $(GCCDB_OBJECTS) -o $(EXECUTABLE)

$(GCCDB_OBJECTS): obj/%_gdb.o: src/%.cpp
	$(GCC) -c -o $@ $(GCC_DBCFLAGS) $<

intel-debug: $(ICCDB_OBJECTS)
	$(ICC) $(ICC_DBLFLAGS) $(LIBS) $(ICCDB_OBJECTS) -o $(EXECUTABLE)

$(ICCDB_OBJECTS): obj/%_idb.o: src/%.cpp
	$(ICC) -c -o $@ $(ICC_DBCFLAGS) $<

pathscale-debug: $(ICCDB_OBJECTS)
	$(PCC) $(PCC_DBLFLAGS) $(LIBS) $(PCCDB_OBJECTS) -o $(EXECUTABLE)

$(PCCDB_OBJECTS): obj/%_pdb.o: src/%.cpp
	$(PCC) -c -o $@ $(PCC_DBCFLAGS) $<

#ibm-debug: $(ICCDB_OBJECTS)
#        $(PCC) $(PCC_DBLFLAGS) $(PCCDB_OBJECTS) -o $(EXECUTABLE)

#$(IBMDB_OBJECTS): obj/%_pdb.o: src/%.cpp
#	$(PCC) -c -o $@ $(PCC_DBCFLAGS) $<

# MPI Targets

mpi-gcc: $(MPI_OBJECTS)
	$(MPICC) $(GCC_LDFLAGS) $(LIBS) $(MPI_OBJECTS) -o $(EXECUTABLE)
#export OMPI_CXX=icc
$(MPI_OBJECTS): obj/%_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(GCC_CFLAGS) $<

mpi-intel: $(MPI_ICC_OBJECTS)
	$(MPICC) $(ICC_LDFLAGS) $(LIBS) $(MPI_ICC_OBJECTS) -o $(EXECUTABLE)
$(MPI_ICC_OBJECTS): obj/%_i_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(ICC_CFLAGS) $<

mpi-pathscale: $(MPI_PCC_OBJECTS)
	$(MPICC) $(PCC_LDFLAGS) $(LIBS) $(MPI_PCC_OBJECTS) -o $(EXECUTABLE)
$(MPI_PCC_OBJECTS): obj/%_p_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(PCC_CFLAGS) $<

mpi-ibm: $(MPI_IBM_OBJECTS)
	$(MPICC) $(IBM_LDFLAGS) $(MPI_IBM_OBJECTS) -o $(EXECUTABLE)
$(MPI_IBM_OBJECTS): obj/%_ibm_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(IBM_CFLAGS) $<

mpi-gcc-debug: $(MPI_GCCDB_OBJECTS)
	$(MPICC) $(GCC_DBLFLAGS) $(LIBS) $(MPI_GCCDB_OBJECTS) -o $(EXECUTABLE)

$(MPI_GCCDB_OBJECTS): obj/%_gdb_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(GCC_DBCFLAGS) $<

mpi-intel-debug: $(MPI_ICCDB_OBJECTS)
	$(MPICC) $(ICC_DBLFLAGS) $(LIBS) $(MPI_ICCDB_OBJECTS) -o $(EXECUTABLE)

$(MPI_ICCDB_OBJECTS): obj/%_idb_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(ICC_DBCFLAGS) $<

mpi-pathscale-debug: $(MPI_PCCDB_OBJECTS)
	$(MPICC) $(PCC_DBLFLAGS) $(LIBS) $(MPI_PCCDB_OBJECTS) -o $(EXECUTABLE)

$(MPI_PCCDB_OBJECTS): obj/%_pdb_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(PCC_DBCFLAGS) $<

# cuda targets
gcc-cuda: obj/cuda/LLG_cuda.o $(CUDA_OBJECTS) 
	$(ICC) $(ICC_LDFLAGS) $(LIBS)  $(CUDALIBS) $(CUDA_OBJECTS) obj/cuda/LLG_cuda.o -o $(EXECUTABLE)

$(CUDA_OBJECTS): obj/%_cuda.o: src/%.cpp
	$(ICC) -c -o $@ $(ICC_CFLAGS) -DCUDA $<

obj/cuda/LLG_cuda.o : src/cuda/LLG_cuda.cu
	nvcc -I/usr/local/cuda/include -I./hdr --compiler-bindir=/usr/bin/g++-4.2 --compiler-options=-O3,-DCUDA  --ptxas-options=-v --maxrregcount=32 -arch=sm_13 -O3  -c $< -o $@

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

package:
	@bash .pack


