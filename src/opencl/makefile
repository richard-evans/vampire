#--------------------------------------------------------------
#          Makefile for vopencl module
#--------------------------------------------------------------

LD=ld

OPENCL_CFLAGS=-Wall -Wextra -std=c++11 -pedantic -Wno-ignored-attributes -Wno-unused-function -DOPENCL_DEBUG

OPENCL_LIBS=-lOpenCL
OPENCL_PR=

# List module object filenames
OPENCL_OBJECTS =\
obj/opencl/config.o \
obj/opencl/data.o \
obj/opencl/dipole.o \
obj/opencl/exchange_fields.o \
obj/opencl/external_fields.o \
obj/opencl/finalize.o \
obj/opencl/initialize.o \
obj/opencl/llg_heun.o \
obj/opencl/opencl_utils.o \
obj/opencl/random.o \
obj/opencl/spin_fields.o \
obj/opencl/statistics.o

OPENCL_OBJ=obj/opencl/vampire_opencl.o

$(OPENCL_OBJ): $(OPENCL_OBJECTS)
	$(LD) -r $^ -o $@

$(OPENCL_OBJECTS): obj/%.o: src/%.cpp
	$(GCC) $< -c -g -o $@ $(GCC_CFLAGS) $(OPENCL_PR) $(OPENCL_LIBS) $(OPENCL_CFLAGS)

gcc-opencl: GCC_CFLAGS+=-DOPENCL
gcc-opencl: OPENCL_PR+=-DOPENCL_DP
gcc-opencl: $(OBJECTS) $(OPENCL_OBJ)
	$(GCC) $^ $(GCC_LDFLAGS) $(LIBS) $(OPENCL_LIBS) -o $(EXECUTABLE)

gcc-opencl-sp: GCC_CFLAGS+=-DOPENCL
gcc-opencl-sp: OPENCL_PR+=-DOPENCL_SP
gcc-opencl-sp: $(OBJECTS) $(OPENCL_OBJ)
	$(GCC) $^ $(GCC_LDFLAGS) $(LIBS) $(OPENCL_LIBS) -o $(EXECUTABLE)

gcc-opencl-native: GCC_CFLAGS+=-DOPENCL_USE_NATIVE_FUNCTIONS
gcc-opencl-native: gcc-opencl

gcc-opencl-native-sp: GCC_CFLAGS+=-DOPENCL_USE_NATIVE_FUNCTIONS
gcc-opencl-native-sp: gcc-opencl-sp
