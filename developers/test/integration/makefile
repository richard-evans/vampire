#===================================================================
#
#             Makefile for Vampire Test Suite
#
#===================================================================

# Compilers
GCC=g++

# LIBS
LIBS= #-lstdc++

# Flags
GCC_CFLAGS=-O3 -std=c++17

# Objects
OBJECTS= \
obj/main.o \
obj/checkpoint.o \
obj/convergence.o \
obj/energy.o \
obj/integrator.o \
obj/montecarlo.o \
obj/structure.o \
obj/utilities.o

EXECUTABLE=integration_tests

all: $(OBJECTS) gcc

# Serial Targets
gcc: $(OBJECTS)
	$(GCC) $(OBJECTS) $(GCC_CFLAGS) $(LIBS) -o $(EXECUTABLE)

$(OBJECTS): obj/%.o: ./src/%.cpp
	@mkdir -p $(@D)
	$(GCC) -c -o $@ $(GCC_CFLAGS) $<

clean:
	@rm -f $(OBJECTS)

purge: clean
	@rm -f $(EXECUTABLE)
