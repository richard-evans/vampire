#--------------------------------------------------------------
#          Makefile for vopencl module
#--------------------------------------------------------------

# Compile with `make opencl(-debug)` or `OPENCL_CFLAGS="options" make opencl(-debug)`
# where available options are:
#
# -DENABLE_MULTIPLE_DEVICES:
#      Enable the ability to use two devices to execute kernels. By
#      default only one is used. To use two devices at least one of
#      gpu:platform-other=n or gpu:device-other=n must be specified in the
#      input file. This will add extra overhead even if only one device is
#      used and presently gives worse performance most of the time even on
#      two devices.
#
# -DOPENCL_TIME_KERNELS:
#      Enable timing of kernels to measure individual kernel performance. This
#      introduces more synchronisation and therefore will result in lower
#      performance overall.
#
# -DOPENCL_USE_NATIVE_FUNCTIONS:
#      Use native device functions which may give better performance but at the
#      cost of less accuracy.
#
# -DOPENCL_USE_VECTOR_TYPE:
#      OpenCL 3-component vector types will be used for data storage. This may
#      (depending on implementation and hardware) increase amount of device
#      memory used but memory accesses may be significantly faster.
#
# -UOPENCL_DP:
#      Use single precision floating point math instead of double precision. It
#      should give better performance and some devices may require this.
#


OPENCL_CFLAGS:=\
-MMD \
-Wno-ignored-attributes -Wno-unused-function \
-DNDEBUG -DOPENCL_DP -DOPENCL_LOG \
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
