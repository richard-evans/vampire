#--------------------------------------------------------------
#                 Makefile for cuda module
#--------------------------------------------------------------




# List module object filenames
cuda_objects= \
	     obj/cuda/data.o \
	     obj/cuda/internal.o \
	     obj/cuda/initialize.o \
	     obj/cuda/fields.o \
	     obj/cuda/llg_heun.o \
	     obj/cuda/finalize.o \
	     obj/cuda/statistics.o

# The device linked object for CUDA
CUDA_DLINK= obj/cuda/vampire_cuda.o

CUDA_CFLAGS= -I./hdr -DCUDA -rdc=true

$(cuda_objects): obj/%.o: src/%.cu
	nvcc -arch=sm_20 $< -c -o $@ -DCUDA $(CUDA_CFLAGS)

$(CUDA_DLINK): $(cuda_objects)
	nvcc -arch=sm_20 -dlink $(cuda_objects) -o $@

# Target to link the CUDA enabled vampire
gcc-cuda: GCC_CFLAGS += -DCUDA
gcc-cuda: $(OBJECTS) $(cuda_objects)
	nvcc $(GCC_LDFLAGS) $(LIBS) ${CUDALIBS} $(OBJECTS) $(cuda_objects) -o $(EXECUTABLE)

# CUDA debug target
gcc-cuda-debug: CUDA_CFLAGS += -DCUDA_DEBUG
gcc-cuda-debug: gcc-cuda
