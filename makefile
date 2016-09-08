#===================================================================
#
#								Makefile for VAMPIRE
#
#===================================================================

export OMPI_CXX=CC
#export OMPI_CXX=icc
#export OMPI_CXX=pathCC
#export MPICH_CXX=g++
export MPICH_CXX=bgxlc++
# Compilers
ICC=icc -DCOMP='"Intel C++ Compiler"'
GCC=g++ -DCOMP='"GNU C++ Compiler"'
LLVM=g++ -DCOMP='"LLVM C++ Compiler"'
PCC=pathCC -DCOMP='"Pathscale C++ Compiler"'
IBM=bgxlc++ -DCOMP='"IBM XLC++ Compiler"'
#MPICC=mpicxx -DMPICF
MPICC=CC -DMPICF

CCC_CFLAGS=-I./hdr -I./src/qvoronoi -O0
CCC_LDFLAGS=-I./hdr -I./src/qvoronoi -O0

export LANG=C
export LC_ALL=C

# LIBS
LIBS=
#-lstdc++
CUDALIBS=-L/usr/local/cuda/lib64/ -lcuda -lcudart
# Debug Flags
ICC_DBCFLAGS= -O0 -C -I./hdr -I./src/qvoronoi
ICC_DBLFLAGS= -C -I./hdr -I./src/qvoronoi

GCC_DBCFLAGS= -Wall -Wextra -O0 -fbounds-check -pedantic -std=c++98 -Wno-long-long -I./hdr -I./src/qvoronoi
GCC_DBLFLAGS= -lstdc++ -fbounds-check -I./hdr -I./src/qvoronoi

PCC_DBCFLAGS= -O0 -I./hdr -I./src/qvoronoi
PCC_DBLFLAGS= -O0 -I./hdr -I./src/qvoronoi

IBM_DBCFLAGS= -O0 -Wall -pedantic -Wextra -I./hdr -I./src/qvoronoi
IBM_DBLFLAGS= -O0 -Wall -pedantic -Wextra -I./hdr -I./src/qvoronoi

# Performance Flags
ICC_CFLAGS= -O3 -axSSE3 -fno-alias -align -falign-functions -I./hdr -I./src/qvoronoi
ICC_LDFLAGS= -I./hdr -I./src/qvoronoi -axSSE3
#ICC_CFLAGS= -O3 -xT -ipo -static -fno-alias -align -falign-functions -vec-report -I./hdr
#ICC_LDFLAGS= -lstdc++ -ipo -I./hdr -xT -vec-report

LLVM_CFLAGS= -O3 -mtune=native -funroll-loops -I./hdr -I./src/qvoronoi
LLVM_LDFLAGS= -lstdc++ -I./hdr -I./src/qvoronoi

GCC_CFLAGS=-O3 -mtune=native -funroll-all-loops -fexpensive-optimizations -funroll-loops -I./hdr -I./src/qvoronoi
GCC_LDFLAGS= -lstdc++ -I./hdr -I./src/qvoronoi

PCC_CFLAGS=-O2 -march=barcelona -ipa -I./hdr -I./src/qvoronoi
PCC_LDFLAGS= -I./hdr -I./src/qvoronoi -O2 -march=barcelona -ipa

NVCC_FLAGS=-I/usr/local/cuda/include -I./hdr -I./src/qvoronoi --compiler-bindir=/usr/bin/g++-4.2 --compiler-options=-O3,-DCUDA
--ptxas-options=-v --maxrregcount=32 -arch=sm_13 -O3
NVCC=nvcc -DCOMP='"GNU C++ Compiler"'

IBM_CFLAGS=-O5 -qarch=450 -qtune=450 -I./hdr -I./src/qvoronoi
IBM_LDFLAGS= -lstdc++ -I./hdr -I./src/qvoronoi -O5 -qarch=450 -qtune=450

CRAY_CFLAGS= -O3 -hfp3 -I./hdr -I./src/qvoronoi
CRAY_LDFLAGS= -I./hdr -I./src/qvoronoi

# Objects
OBJECTS= \
obj/create/create_system2.o \
obj/create/cs_create_crystal_structure2.o \
obj/create/cs_create_system_type2.o \
obj/create/cs_create_neighbour_list2.o \
obj/create/cs_particle_shapes.o \
obj/create/cs_set_atom_vars2.o \
obj/create/cs_voronoi2.o \
obj/create/multilayers.o \
obj/data/atoms.o \
obj/data/category.o \
obj/data/cells.o \
obj/data/grains.o \
obj/data/lattice_anisotropy.o \
obj/main/initialise_variables.o \
obj/main/main.o \
obj/main/material.o \
obj/mpi/LLGHeun-mpi.o \
obj/mpi/LLGMidpoint-mpi.o \
obj/mpi/mpi_generic.o \
obj/mpi/mpi_create2.o \
obj/mpi/mpi_comms.o \
obj/mpi/wrapper.o \
obj/program/bmark.o \
obj/program/cmc_anisotropy.o \
obj/program/curie_temperature.o \
obj/program/diagnostics.o \
obj/program/field_cool.o \
obj/program/hamr.o \
obj/program/hybrid_cmc.o \
obj/program/hysteresis.o \
obj/program/lagrange.o \
obj/program/LLB_Boltzmann.o \
obj/program/partial_hysteresis.o \
obj/program/static_hysteresis.o \
obj/program/time_series.o \
obj/program/temperature_pulse.o \
obj/program/localised_temperature_pulse.o \
obj/program/effective_damping.o \
obj/program/fmr.o \
obj/random/mtrand.o \
obj/random/random.o \
obj/simulate/energy.o \
obj/simulate/fields.o \
obj/simulate/demag.o \
obj/simulate/LLB.o \
obj/simulate/LLGHeun.o \
obj/simulate/LLGMidpoint.o \
obj/simulate/mc.o \
obj/simulate/mc_moves.o \
obj/simulate/cmc.o \
obj/simulate/cmc_mc.o \
obj/simulate/sim.o \
obj/simulate/standard_programs.o \
obj/statistics/data.o \
obj/statistics/initialize.o \
obj/statistics/magnetization.o \
obj/statistics/statistics.o \
obj/statistics/susceptibility.o \
obj/utility/checkpoint.o \
obj/utility/errors.o \
obj/utility/statistics.o \
obj/utility/units.o \
obj/utility/vconfig.o \
obj/utility/vio.o \
obj/utility/vmath.o\
obj/qvoronoi/geom.o\
obj/qvoronoi/geom2.o\
obj/qvoronoi/global.o\
obj/qvoronoi/io.o\
obj/qvoronoi/libqhull.o\
obj/qvoronoi/mem.o\
obj/qvoronoi/merge.o\
obj/qvoronoi/poly.o\
obj/qvoronoi/poly2.o\
obj/qvoronoi/qhrandom.o\
obj/qvoronoi/qset.o\
obj/qvoronoi/qvoronoi.o\
obj/qvoronoi/rboxlib.o\
obj/qvoronoi/stat.o\
obj/qvoronoi/user.o\
obj/qvoronoi/usermem.o\
obj/qvoronoi/userprintf.o\
obj/qvoronoi/userprintf_rbox.o\

# Include supplementary makefiles
include src/create/makefile
include src/gpu/makefile
include src/ltmp/makefile
include src/simulate/makefile

ICC_OBJECTS=$(OBJECTS:.o=_i.o)
LLVM_OBJECTS=$(OBJECTS:.o=_llvm.o)
IBM_OBJECTS=$(OBJECTS:.o=_ibm.o)
ICCDB_OBJECTS=$(OBJECTS:.o=_idb.o)
GCCDB_OBJECTS=$(OBJECTS:.o=_gdb.o)
PCCDB_OBJECTS=$(OBJECTS:.o=_pdb.o)
IBMDB_OBJECTS=$(OBJECTS:.o=_ibmdb.o)

MPI_OBJECTS=$(OBJECTS:.o=_mpi.o)
MPI_ICC_OBJECTS=$(OBJECTS:.o=_i_mpi.o)
MPI_LLVM_OBJECTS=$(OBJECTS:.o=_llvm_mpi.o)
MPI_PCC_OBJECTS=$(OBJECTS:.o=_p_mpi.o)
MPI_IBM_OBJECTS=$(OBJECTS:.o=_ibm_mpi.o)
MPI_CRAY_OBJECTS=$(OBJECTS:.o=_cray_mpi.o)
MPI_ICCDB_OBJECTS=$(OBJECTS:.o=_idb_mpi.o)
MPI_GCCDB_OBJECTS=$(OBJECTS:.o=_gdb_mpi.o)
MPI_PCCDB_OBJECTS=$(OBJECTS:.o=_pdb_mpi.o)
MPI_IBMDB_OBJECTS=$(OBJECTS:.o=_ibmdb_mpi.o)
MPI_CRAY_OBJECTS=$(OBJECTS:.o=_cray_mpi.o)

CUDA_OBJECTS=$(OBJECTS:.o=_cuda.o)
EXECUTABLE=vampire

all: $(OBJECTS) serial

# Serial Targets
serial: $(OBJECTS)
	$(GCC) $(GCC_LDFLAGS) $(LIBS) $(OBJECTS) -o $(EXECUTABLE)

$(OBJECTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $<

serial-intel: $(ICC_OBJECTS)
	$(ICC) $(ICC_LDFLAGS) $(LIBS) $(ICC_OBJECTS) -o $(EXECUTABLE)

$(ICC_OBJECTS): obj/%_i.o: src/%.cpp
	$(ICC) -c -o $@ $(ICC_CFLAGS) $<

serial-llvm: $(LLVM_OBJECTS)
	$(LLVM) $(LLVM_LDFLAGS) $(LIBS) $(LLVM_OBJECTS) -o $(EXECUTABLE)

$(LLVM_OBJECTS): obj/%_llvm.o: src/%.cpp
	$(LLVM) -c -o $@ $(LLVM_CFLAGS) $<

serial-ibm: $(IBM_OBJECTS)
	$(IBM) $(IBM_LDFLAGS) $(IBM_OBJECTS) -o $(EXECUTABLE)

$(IBM_OBJECTS): obj/%_ibm.o: src/%.cpp
	$(IBM) -c -o $@ $(IBM_CFLAGS) $<

serial-debug: $(GCCDB_OBJECTS)
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

parallel: $(MPI_OBJECTS)
	$(MPICC) $(GCC_LDFLAGS) $(LIBS) $(MPI_OBJECTS) -o $(EXECUTABLE)
#export OMPI_CXX=icc
$(MPI_OBJECTS): obj/%_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(GCC_CFLAGS) $<

parallel-intel: $(MPI_ICC_OBJECTS)
	$(MPICC) $(ICC_LDFLAGS) $(LIBS) $(MPI_ICC_OBJECTS) -o $(EXECUTABLE)
$(MPI_ICC_OBJECTS): obj/%_i_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(ICC_CFLAGS) $<intel: $(MPI_ICC_OBJECTS)

parallel-cray: $(MPI_CRAY_OBJECTS)
	$(MPICC) $(CRAY_LDFLAGS) $(LIBS) $(MPI_CRAY_OBJECTS) -o $(EXECUTABLE)
$(MPI_CRAY_OBJECTS): obj/%_cray_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(CRAY_CFLAGS) $<

parallel-llvm: $(MPI_LLVM_OBJECTS)
	$(MPICC) $(LLVM_LDFLAGS) $(LIBS) $(MPI_LLVM_OBJECTS) -o $(EXECUTABLE)
$(MPI_LLVM_OBJECTS): obj/%_llvm_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(LLVM_CFLAGS) $<

parallel-pathscale: $(MPI_PCC_OBJECTS)
	$(MPICC) $(PCC_LDFLAGS) $(LIBS) $(MPI_PCC_OBJECTS) -o $(EXECUTABLE)
$(MPI_PCC_OBJECTS): obj/%_p_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(PCC_CFLAGS) $<

parallel-ibm: $(MPI_IBM_OBJECTS)
	$(MPICC) $(IBM_LDFLAGS) $(MPI_IBM_OBJECTS) -o $(EXECUTABLE)
$(MPI_IBM_OBJECTS): obj/%_ibm_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(IBM_CFLAGS) $<

parallel-debug: $(MPI_GCCDB_OBJECTS)
	$(MPICC) $(GCC_DBLFLAGS) $(LIBS) $(MPI_GCCDB_OBJECTS) -o $(EXECUTABLE)

$(MPI_GCCDB_OBJECTS): obj/%_gdb_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(GCC_DBCFLAGS) $<

parallel-intel-debug: $(MPI_ICCDB_OBJECTS)
	$(MPICC) $(ICC_DBLFLAGS) $(LIBS) $(MPI_ICCDB_OBJECTS) -o $(EXECUTABLE)

$(MPI_ICCDB_OBJECTS): obj/%_idb_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(ICC_DBCFLAGS) $<

parallel-cray: $(MPI_CRAY_OBJECTS)
	$(MPICC) $(CCC_LDFLAGS) $(LIBS) $(MPI_CRAY_OBJECTS) -o $(EXECUTABLE)
$(MPI_CRAY_OBJECTS): obj/%_cray_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(CCC_CFLAGS) $<

parallel-pathscale-debug: $(MPI_PCCDB_OBJECTS)
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
	@rm -f hdr/*~
	@rm -f src/*~
	@rm -f src/*/*~
