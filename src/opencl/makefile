#--------------------------------------------------------------
#          Makefile for vopencl module
#--------------------------------------------------------------

# Compile as `make opencl(-debug)` or `OPENCL_CFLAGS="options" make opencl(-debug)`
# where available options are:
#
# -DOPENCL_USE_NATIVE_FUNCTIONS:
#      Use native device functions which may give better performance but at the cost
#      of less accuracy.
#
# -UOPENCL_DP:
#      Use single precision floating point math instead of double precision. Some
#      devices may require this.
#
# -UUSE_VECTOR_TYPE:
#      OpenCL vector types will not be used for data storage. This may (depending on implementation
#      and hardware) decrease amount of device memory used but memory accesses may be significantly
#      slower.


OPENCL_CFLAGS:=\
-MMD \
-Wno-ignored-attributes -Wno-unused-function \
-DNDEBUG -DOPENCL_DP -DUSE_VECTOR_TYPE -DOPENCL_LOG \
$(OPENCL_CFLAGS)

OPENCL_DEBUG_CFLAGS:=\
-Og -g \
-fno-fast-math -Wall -Wextra -Wpedantic \
-UNDEBUG -DOPENCL_DEBUG -DOPENCL_TIME_KERNELS

OPENCL_LIBS=-lOpenCL

# List module object filenames
OPENCL_OBJECTS =\
obj/opencl/build_kernels.o \
obj/opencl/config.o \
obj/opencl/data.o \
obj/opencl/dipole.o \
obj/opencl/finalize.o \
obj/opencl/initialize.o \
obj/opencl/init_kernels.o \
obj/opencl/llg_heun.o \
obj/opencl/opencl_utils.o \
obj/opencl/statistics.o \

# each object depends on .cpp and .hpp
# this will rebuild them if the .hpp changes
# generated with -MMD cflag
OPENCL_DEPS=$(OPENCL_OBJECTS:%.o=%.d)
-include $(OPENCL_DEPS)

# target to delete only the OpenCL module objects and dependencies
opencl-clean:
	rm -f $(OPENCL_OBJECTS) $(OPENCL_DEPS)

opencl: GCC_CFLAGS+=-std=c++11 -DOPENCL
opencl: $(OBJECTS) $(OPENCL_OBJECTS)
	$(GCC) $^ $(GCC_LDFLAGS) $(LIBS) $(OPENCL_LIBS) -o $(EXECUTABLE)

opencl-debug:
	OPENCL_CFLAGS+=" $(OPENCL_DEBUG_CFLAGS)" $(MAKE) opencl



obj/opencl/%.o: src/opencl/%.cpp
	$(GCC) $< -c -o $@ $(GCC_CFLAGS) $(OPENCL_CFLAGS)
