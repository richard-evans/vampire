#--------------------------------------------------------------
#          Makefile for vopencl module
#--------------------------------------------------------------

LD=ld

OPENCL_LIBS=-lOpenCL
OPENCL_PR=-DOPENCL_DP

# List module object filenames
OPENCL_OBJECTS =\
obj/opencl/config.o \
obj/opencl/data.o \
obj/opencl/dipole.o \
obj/opencl/exchange_fields.o \
obj/opencl/external_fields.o \
obj/opencl/finalize.o \
obj/opencl/initialize.o \
obj/opencl/interface.o \
obj/opencl/llg_heun.o \
obj/opencl/opencl_utils.o \
obj/opencl/random.o \
obj/opencl/spin_fields.o \
obj/opencl/statistics.o

OPENCL_OBJ=obj/opencl/vampire_opencl.o

$(OPENCL_OBJ): $(OPENCL_OBJECTS)
	$(LD) -r $^ -o $@

$(OPENCL_OBJECTS): obj/%.o: src/%.cpp
	$(GCC) $< -c -g -o $@ $(GCC_CFLAGS) $(OPENCL_PR) $(OPENCL_LIBS) -DOPENCL_DEBUG -Wall -Wextra -pedantic

gcc-opencl: GCC_CFLAGS+=-DOPENCL -std=c++11
gcc-opencl: $(OBJECTS) $(OPENCL_OBJ)
	$(GCC) $^ $(GCC_LDFLAGS) $(LIBS) $(OPENCL_LIBS) -o $(EXECUTABLE)
