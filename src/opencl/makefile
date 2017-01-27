#--------------------------------------------------------------
#          Makefile for vopencl module
#--------------------------------------------------------------

OPENCL_CFLAGS=\
-std=c++11 -pedantic \
-Wall -Wextra -Wno-ignored-attributes -Wno-unused-function \
-MMD -DOPENCL_DEBUG -DUSE_VECTOR_TYPE

OPENCL_LIBS=-lOpenCL
OPENCL_PR=

# List module object filenames
OPENCL_OBJECTS =\
obj/opencl/config.o \
obj/opencl/data.o \
obj/opencl/dipole.o \
obj/opencl/finalize.o \
obj/opencl/initialize.o \
obj/opencl/init_kernels.o \
obj/opencl/llg_heun.o \
obj/opencl/opencl_utils.o \
obj/opencl/statistics.o

# each object depends on .cpp and .hpp
# this will rebuild them if the .hpp changes
# generated with -MMD cflag
-include $(OPENCL_OBJECTS:%.o=%.d)

OPENCL_OBJ=obj/opencl/vampire_opencl.o

# double precision
gcc-opencl: GCC_CFLAGS+=-DOPENCL
gcc-opencl: OPENCL_PR+=-DOPENCL_DP
gcc-opencl: $(OBJECTS) $(OPENCL_OBJ)
	$(GCC) $^ $(GCC_LDFLAGS) $(LIBS) $(OPENCL_LIBS) -o $(EXECUTABLE)

# single precision
gcc-opencl-sp: GCC_CFLAGS+=-DOPENCL
gcc-opencl-sp: OPENCL_PR+=-DOPENCL_SP
gcc-opencl-sp: $(OBJECTS) $(OPENCL_OBJ)
	$(GCC) $^ $(GCC_LDFLAGS) $(LIBS) $(OPENCL_LIBS) -o $(EXECUTABLE)

# double precision with native functions
gcc-opencl-native: GCC_CFLAGS+=-DOPENCL_USE_NATIVE_FUNCTIONS
gcc-opencl-native: gcc-opencl

# single precision with native functions
gcc-opencl-native-sp: GCC_CFLAGS+=-DOPENCL_USE_NATIVE_FUNCTIONS
gcc-opencl-native-sp: gcc-opencl-sp


$(OPENCL_OBJ): $(OPENCL_OBJECTS)
	$(LD) -r $^ -o $@

obj/opencl/%.o: src/opencl/%.cpp
	$(GCC) $< -c -g -o $@ $(GCC_CFLAGS) $(OPENCL_PR) $(OPENCL_LIBS) $(OPENCL_CFLAGS)
