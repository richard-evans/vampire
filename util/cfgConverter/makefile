#===================================================================
#
#Makefile for cfgConverter
#
#===================================================================

# Compilers
GCC=g++ -DCOMP='"GNU C++ Compiler"'

# LIBS
LIBS=-lstdc++ 

# Flags
GCC_CFLAGS=-O3 -mtune=native -funroll-all-loops -fexpensive-optimizations -funroll-loops -std=c++11 

# Objects
OBJECTS= \
../../obj/util/colour.o \
../../obj/util/manipulate.o \
../../obj/util/main.o \
../../obj/util/read.o \
../../obj/util/output.o \


EXECUTABLE=cfgConverter

all: $(OBJECTS) gcc

# Serial Targets
gcc: $(OBJECTS)
	$(GCC) $(OBJECTS) $(LIBS) -o $(EXECUTABLE)

$(OBJECTS): ../../obj/util/%.o: ./%.cpp
	$(GCC) -c -o $@ $(GCC_CFLAGS) $<

clean:
	@rm -f obj/*.o
	@rm -f obj/*/*.o

purge:
	@rm -f obj/*.o
	@rm -f obj/*/*.o
	@rm -f $(EXECUTABLE)
