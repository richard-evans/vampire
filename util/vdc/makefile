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
obj/main.o \
obj/non_magnetic.o \
obj/povray.o \
obj/spins.o \
obj/ssc.o \
obj/text.o \
obj/vtk.o \
obj/xyz.o

EXECUTABLE=vdc

all: $(OBJECTS) gcc

# Serial Targets
gcc: $(OBJECTS)
	$(GCC) $(OBJECTS) $(GCC_CFLAGS) $(LIBS) -o $(EXECUTABLE)

$(OBJECTS): obj/%.o: ./%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $<

clean:
	@rm -f obj/*.o

purge:
	@rm -f obj/*.o
	@rm -f $(EXECUTABLE)
