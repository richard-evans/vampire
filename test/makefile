#===================================================================
#
#             Makefile for Vampire Test Suite
#
#===================================================================

# Compilers
GCC=g++

# LIBS
LIBS=-lstdc++

# Flags
GCC_CFLAGS=-O3 -std=c++17

# Objects
OBJECTS= \
obj/test_suite.o

EXECUTABLE=tests

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
