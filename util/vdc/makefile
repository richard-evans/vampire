#===================================================================
#
#             Makefile for Vampire Data Converter (vdc)
#
#===================================================================

# Compilers
GCC=g++

# LIBS
LIBS=-lstdc++

# Flags
GCC_CFLAGS=-O3 -std=c++0x
#GCC_CFLAGS=-O3 -std=c++0x -fopenmp -funroll-loops -mavx2

# Objects
OBJECTS= \
obj/cells.o \
obj/colour.o \
obj/colourmaps.o \
obj/command.o \
obj/coords.o \
obj/data.o \
obj/functions.o \
obj/main.o \
obj/non_magnetic.o \
obj/povray.o \
obj/povray-cells.o \
obj/read.o \
obj/spins.o \
obj/ssc.o \
obj/text.o \
obj/vtk.o \
obj/xyz.o

EXECUTABLE=vdc

all: $(OBJECTS) gcc

GCCDB_OBJECTS=$(OBJECTS:.o=_gdb.o)

GCC_DBCFLAGS= -g -Wall -Wextra -O0 -fbounds-check -pedantic -std=c++11 -Wno-long-long
GCC_DBLFLAGS= -g -std=c++11 -fbounds-check

OMP_FLAGS= -fopenmp

# Serial Targets
gcc: $(OBJECTS)
	$(GCC) $(OBJECTS) $(GCC_CFLAGS) $(LIBS) -o $(EXECUTABLE)

gcc-debug: $(GCCDB_OBJECTS)
	$(GCC) $(GCC_DBLFLAGS) $(OMP_FLAGS) $(GCCDB_OBJECTS) -o $(EXECUTABLE)-debug

$(GCCDB_OBJECTS): obj/%_gdb.o: ./%.cpp
	$(GCC) -c -o $@ $(GCC_DBCFLAGS) $(OMP_FLAGS) $<

$(OBJECTS): obj/%.o: ./%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $<

clean:
	@rm -f obj/*.o

purge:
	@rm -f obj/*.o
	@rm -f $(EXECUTABLE) $(EXECUTABLE)-debug
