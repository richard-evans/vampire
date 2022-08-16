#===================================================================
#
#								Makefile for VAMPIRE
#
#===================================================================

# Specify compiler for MPI compilation with openmpi
export OMPI_CXX=g++ -std=c++11

#export OMPI_CXX=icc
#export OMPI_CXX=pathCC
# Specify compiler for MPI compilation with mpich
#export MPICH_CXX=g++
#export MPICH_CXX=bgxlc++

# Include the FFTW library by uncommenting the -DFFT (off by default)
#export incFFT= -DFFT -DFFTW_OMP -fopenmp
#export FFTLIBS= -lfftw3_omp -lfftw3

# Compilers
ICC=icc -std=c++11 -DCOMP='"Intel C++ Compiler"'
GCC=g++ -std=c++11 -DCOMP='"GNU C++ Compiler"'
LLVM=g++ -std=c++11 -DCOMP='"LLVM C++ Compiler"'
PCC=pathCC -DCOMP='"Pathscale C++ Compiler"'
IBM=bgxlc++ -DCOMP='"IBM XLC++ Compiler"'
MPICC=mpicxx -DMPICF
MPIICC=mpiicpc -DMPICF

LIBS= -lstdc++
#-lm $(FFTLIBS) -L/opt/local/lib/

CCC_CFLAGS=-I./hdr -I./src/qvoronoi -O0
CCC_LDFLAGS=-I./hdr -I./src/qvoronoi -O0

export LANG=C
export LC_ALL=C

# LIBS

CUDALIBS=-L/usr/local/cuda/lib64/ -lcuda -lcudart

# Debug Flags
ICC_DBCFLAGS= -O0 -C -I./hdr -I./src/qvoronoi
ICC_DBLFLAGS= -C -I./hdr -I./src/qvoronoi

GCC_DBCFLAGS= -g -pg -fprofile-arcs -ftest-coverage -Wall -Wextra -O0 -fbounds-check -pedantic -std=c++0x -Wno-long-long -I./hdr -I./src/qvoronoi -Wsign-compare
GCC_DBLFLAGS= -g -pg -fprofile-arcs -ftest-coverage -lstdc++ -std=c++0x -fbounds-check -I./hdr -I./src/qvoronoi -Wsign-compare

PCC_DBCFLAGS= -O0 -I./hdr -I./src/qvoronoi
PCC_DBLFLAGS= -O0 -I./hdr -I./src/qvoronoi
IBM_DBCFLAGS= -O0 -Wall -pedantic -Wextra -I./hdr -I./src/qvoronoi
IBM_DBLFLAGS= -O0 -Wall -pedantic -Wextra -I./hdr -I./src/qvoronoi

LLVM_DBCFLAGS= -Wall -Wextra -O0 -pedantic -std=c++11 -Wno-long-long -I./hdr -I./src/qvoronoi -Wsign-compare
LLVM_DBLFLAGS= -Wall -Wextra -O0 -lstdc++ -I./hdr -I./src/qvoronoi -Wsign-compare

# Performance Flags
ICC_CFLAGS= -O3 -axCORE-AVX2 -fno-alias -align -falign-functions -I./hdr -I./src/qvoronoi
ICC_LDFLAGS= -I./hdr -I./src/qvoronoi -axCORE-AVX2
#ICC_CFLAGS= -O3 -xT -ipo -static -fno-alias -align -falign-functions -vec-report -I./hdr
#ICC_LDFLAGS= -lstdc++ -ipo -I./hdr -xT -vec-report

LLVM_CFLAGS= -Wall -pedantic -O3 -mtune=native -funroll-loops -I./hdr -I./src/qvoronoi
LLVM_LDFLAGS= -lstdc++ -I./hdr -I./src/qvoronoi

GCC_CFLAGS=-O3 -mtune=native -funroll-all-loops -fexpensive-optimizations -funroll-loops -I./hdr -I./src/qvoronoi -std=c++11 -Wsign-compare
GCC_LDFLAGS= -lstdc++ -I./hdr -I./src/qvoronoi -Wsign-compare

PCC_CFLAGS=-O2 -march=barcelona -ipa -I./hdr -I./src/qvoronoi
PCC_LDFLAGS= -I./hdr -I./src/qvoronoi -O2 -march=barcelona -ipa


IBM_CFLAGS=-O5 -qarch=450 -qtune=450 -I./hdr -I./src/qvoronoi
IBM_LDFLAGS= -lstdc++ -I./hdr -I./src/qvoronoi -O5 -qarch=450 -qtune=450

CRAY_CFLAGS= -O3 -hfp3 -I./hdr -I./src/qvoronoi
CRAY_LDFLAGS= -I./hdr -I./src/qvoronoi


# Save git commit in simple function
GHASH:=$(shell git rev-parse HEAD)
# special options for certain files

OPTIONS=

# Objects
OBJECTS= \
obj/data/atoms.o \
obj/data/category.o \
obj/data/grains.o \
obj/random/mtrand.o \
obj/random/random.o \
obj/simulate/energy.o \
obj/simulate/fields.o \
obj/simulate/LLB.o \
obj/simulate/LLGHeun.o \
obj/simulate/LLGMidpoint.o \
obj/simulate/sim.o \
obj/simulate/standard_programs.o \
obj/spintorque/data.o \
obj/spintorque/field.o \
obj/spintorque/initialise.o \
obj/spintorque/interface.o \
obj/spintorque/magnetization.o \
obj/spintorque/matrix.o \
obj/spintorque/output.o \
obj/spintorque/spinaccumulation.o \
obj/utility/checkpoint.o \
obj/utility/errors.o \
obj/utility/statistics.o \
obj/utility/units.o \
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
include src/anisotropy/makefile
include src/cells/makefile
include src/create/makefile
include src/config/makefile
include src/constants/makefile
include src/dipole/makefile
include src/exchange/makefile
include src/gpu/makefile
include src/hamr/makefile
include src/hierarchical/makefile
include src/ltmp/makefile
include src/main/makefile
include src/montecarlo/makefile
include src/micromagnetic/makefile
include src/mpi/makefile
include src/neighbours/makefile
include src/program/makefile
include src/simulate/makefile
include src/spintransport/makefile
include src/statistics/makefile
include src/unitcell/makefile
include src/vio/makefile
include src/environment/makefile

# Cuda must be last for some odd reason
include src/cuda/makefile
include src/opencl/makefile

ICC_OBJECTS=$(OBJECTS:.o=_i.o)
LLVM_OBJECTS=$(OBJECTS:.o=_llvm.o)
IBM_OBJECTS=$(OBJECTS:.o=_ibm.o)
ICCDB_OBJECTS=$(OBJECTS:.o=_idb.o)
GCCDB_OBJECTS=$(OBJECTS:.o=_gdb.o)
PCCDB_OBJECTS=$(OBJECTS:.o=_pdb.o)
IBMDB_OBJECTS=$(OBJECTS:.o=_ibmdb.o)
LLVMDB_OBJECTS=$(OBJECTS:.o=_llvmdb.o)

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
MPI_CRAYDB_OBJECTS=$(OBJECTS:.o=_craydb_mpi.o)
MPI_ARCHER_OBJECTS=$(OBJECTS:.o=_archer_mpi.o)

CLEXECUTABLE=vampire-opencl
CUDAEXECUTABLE=vampire-cuda
EXECUTABLE=vampire-serial
PEXECUTABLE=vampire-parallel

# Set default make target in GNU make > v3.81
.DEFAULT_GOAL := all

# make serial and parallel versions and utilities
all: serial parallel vdc

# Serial Targets
serial: $(OBJECTS)
	$(GCC) $(GCC_LDFLAGS)  $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

$(OBJECTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $(OPTIONS) $<

serial-intel: $(ICC_OBJECTS)
	$(ICC) $(ICC_LDFLAGS) $(LIBS) $(ICC_OBJECTS) -o $(EXECUTABLE)-intel

$(ICC_OBJECTS): obj/%_i.o: src/%.cpp
	$(ICC) -c -o $@ $(ICC_CFLAGS) $(OPTIONS) $<

serial-llvm: $(LLVM_OBJECTS)
	$(LLVM) $(LLVM_LDFLAGS) $(LIBS) $(LLVM_OBJECTS) -o $(EXECUTABLE)

$(LLVM_OBJECTS): obj/%_llvm.o: src/%.cpp
	$(LLVM) -c -o $@ $(LLVM_CFLAGS) $(OPTIONS) $<

serial-ibm: $(IBM_OBJECTS)
	$(IBM) $(IBM_LDFLAGS) $(IBM_OBJECTS) -o $(EXECUTABLE)

$(IBM_OBJECTS): obj/%_ibm.o: src/%.cpp
	$(IBM) -c -o $@ $(IBM_CFLAGS) $(OPTIONS) $<

serial-debug: $(GCCDB_OBJECTS)
	$(GCC) $(GCC_DBLFLAGS) $(LIBS) $(GCCDB_OBJECTS) -o $(EXECUTABLE)-debug

$(GCCDB_OBJECTS): obj/%_gdb.o: src/%.cpp
	$(GCC) -c -o $@ $(GCC_DBCFLAGS) $(OPTIONS) $<

serial-llvm-debug: $(LLVMDB_OBJECTS)
	$(LLVM) $(LLVM_DBLFLAGS) $(LIBS) $(LLVMDB_OBJECTS) -o $(EXECUTABLE)

$(LLVMDB_OBJECTS): obj/%_llvmdb.o: src/%.cpp
	$(LLVM) -c -o $@ $(LLVM_DBCFLAGS) $(OPTIONS) $<

intel-debug: $(ICCDB_OBJECTS)
	$(ICC) $(ICC_DBLFLAGS) $(LIBS) $(ICCDB_OBJECTS) -o $(EXECUTABLE)-intel-debug

$(ICCDB_OBJECTS): obj/%_idb.o: src/%.cpp
	$(ICC) -c -o $@ $(ICC_DBCFLAGS) $(OPTIONS) $<

pathscale-debug: $(ICCDB_OBJECTS)
	$(PCC) $(PCC_DBLFLAGS) $(LIBS) $(PCCDB_OBJECTS) -o $(EXECUTABLE)

$(PCCDB_OBJECTS): obj/%_pdb.o: src/%.cpp
	$(PCC) -c -o $@ $(PCC_DBCFLAGS) $(OPTIONS) $<

#ibm-debug: $(ICCDB_OBJECTS)
#        $(PCC) $(PCC_DBLFLAGS) $(PCCDB_OBJECTS) -o $(EXECUTABLE)

#$(IBMDB_OBJECTS): obj/%_pdb.o: src/%.cpp
#	$(PCC) -c -o $@ $(PCC_DBCFLAGS) $<

# MPI Targets

parallel: $(MPI_OBJECTS)
	$(MPICC) $(GCC_LDFLAGS) $(MPI_OBJECTS) $(LIBS) -o $(PEXECUTABLE)

$(MPI_OBJECTS): obj/%_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(GCC_CFLAGS) $(OPTIONS) $<

parallel-intel: $(MPI_ICC_OBJECTS)
	$(MPIICC) $(ICC_LDFLAGS) $(LIBS) $(MPI_ICC_OBJECTS) -o $(PEXECUTABLE)-intel

$(MPI_ICC_OBJECTS): obj/%_i_mpi.o: src/%.cpp
	$(MPIICC) -c -o $@ $(ICC_CFLAGS) $<

parallel-cray: $(MPI_CRAY_OBJECTS)
	$(MPICC) $(CRAY_LDFLAGS) $(LIBS) $(MPI_CRAY_OBJECTS) -o $(PEXECUTABLE)

$(MPI_CRAY_OBJECTS): obj/%_cray_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(CRAY_CFLAGS) $<

parallel-archer: $(MPI_ARCHER_OBJECTS)
	CC -DMPICF $(GCC_LDFLAGS) $(LIBS) $(MPI_ARCHER_OBJECTS) -o $(PEXECUTABLE)

$(MPI_ARCHER_OBJECTS): obj/%_archer_mpi.o: src/%.cpp
	CC -DMPICF -c -o $@ $(GCC_CFLAGS) $<

parallel-llvm: $(MPI_LLVM_OBJECTS)
	$(MPICC) $(LLVM_LDFLAGS) $(LIBS) $(MPI_LLVM_OBJECTS) -o $(PEXECUTABLE)

$(MPI_LLVM_OBJECTS): obj/%_llvm_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(LLVM_CFLAGS) $(OPTIONS) $<

parallel-pathscale: $(MPI_PCC_OBJECTS)
	$(MPICC) $(PCC_LDFLAGS) $(LIBS) $(MPI_PCC_OBJECTS) -o $(PEXECUTABLE)

$(MPI_PCC_OBJECTS): obj/%_p_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(PCC_CFLAGS) $(OPTIONS) $<

parallel-ibm: $(MPI_IBM_OBJECTS)
	$(MPICC) $(IBM_LDFLAGS) $(MPI_IBM_OBJECTS) -o $(PEXECUTABLE)

$(MPI_IBM_OBJECTS): obj/%_ibm_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(IBM_CFLAGS) $(OPTIONS) $<

parallel-debug: $(MPI_GCCDB_OBJECTS)
	$(MPICC) $(GCC_DBLFLAGS) $(LIBS) $(MPI_GCCDB_OBJECTS) -o $(PEXECUTABLE)-debug

$(MPI_GCCDB_OBJECTS): obj/%_gdb_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(GCC_DBCFLAGS) $(OPTIONS) $<

parallel-intel-debug: $(MPI_ICCDB_OBJECTS)
	$(MPIICC) $(ICC_DBLFLAGS) $(LIBS) $(MPI_ICCDB_OBJECTS) -o $(PEXECUTABLE)-intel-debug

$(MPI_ICCDB_OBJECTS): obj/%_idb_mpi.o: src/%.cpp
	$(MPIICC) -c -o $@ $(ICC_DBCFLAGS) $(OPTIONS) $<

parallel-cray-debug: $(MPI_CRAY_OBJECTS)
	$(MPICC) $(CCC_LDFLAGS) $(LIBS) $(MPI_CRAYDB_OBJECTS) -o $(PEXECUTABLE)

$(MPI_CRAYDB_OBJECTS): obj/%_craydb_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(CCC_CFLAGS) $(OPTIONS) $<

parallel-pathscale-debug: $(MPI_PCCDB_OBJECTS)
	$(MPICC) $(PCC_DBLFLAGS) $(LIBS) $(MPI_PCCDB_OBJECTS) -o $(PEXECUTABLE)

$(MPI_PCCDB_OBJECTS): obj/%_pdb_mpi.o: src/%.cpp
	$(MPICC) -c -o $@ $(PCC_DBCFLAGS) $(OPTIONS) $<

clean:
	@rm -f obj/*.o
	@rm -f obj/*/*.o

purge:
	@rm -f obj/*.o
	@rm -f obj/*/*.o
	@rm -f vampire-*

tidy:
	@rm -f *~
	@rm -f hdr/*~
	@rm -f src/*~
	@rm -f src/*/*~

tests:
	$(MAKE) -C test/integration/
	$(MAKE) -C test/unit/

vdc:
	$(MAKE) -C util/vdc/

vdc-debug:
	$(MAKE) -C util/vdc/ gcc-debug

vdc-purge:
	$(MAKE) -C util/vdc/ purge

install:
	echo "Preparing installation package"
	rm -rf vampire.pkg
	mkdir vampire.pkg
	mkdir vampire.pkg/bin
	cp vampire-* vampire.pkg/bin/
	cp util/vdc/vdc vampire.pkg/bin/
	mkdir vampire.pkg/examples
	cp input vampire.pkg/examples/
	cp Co.mat vampire.pkg/examples/
	sudo mv -f vampire.pkg /opt/vampire
	sudo echo "/opt/vampire/bin" > /etc/paths.d/vampire_path

uninstall:
	rm -rf /opt/vampire
	rm -f /etc/paths.d/vampire_path
