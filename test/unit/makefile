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
GCC_CFLAGS=-O3 -std=c++17 -I../../hdr/

# Objects
TEST_OBJECTS= \
obj/unit_tests.o \
obj/utility/units_test.o \
obj/utility/utility_test.o

VAMPIRE_OBJECTS= \
../../obj/main/githash.o \
../../obj/main/version.o \
../../obj/utility/errors.o \
../../obj/utility/units.o \
../../obj/vio/data.o \
../../obj/vio/globalio.o \
../../obj/vio/timestamp.o

EXECUTABLE=unit_tests

all: $(TEST_OBJECTS) $(VAMPIRE_OBJECTS) gcc

# Serial Targets
gcc: $(TEST_OBJECTS) $(VAMPIRE_OBJECTS)
	$(GCC) $(TEST_OBJECTS) $(VAMPIRE_OBJECTS) $(GCC_CFLAGS) $(LIBS) -o $(EXECUTABLE)

$(TEST_OBJECTS): obj/%.o: src/%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $<

# compile vampire objects to vampire object folder
$(VAMPIRE_OBJECTS): ../../obj/%.o: ../../src/%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $<

clean:
	@rm -f obj/*.o
	@rm -f ../../obj/*.o

purge:
	@rm -f obj/*.o
	@rm -f $(EXECUTABLE)
