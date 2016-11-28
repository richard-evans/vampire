#--------------------------------------------------------------
#                 Makefile for cuda module
#--------------------------------------------------------------

NVCC_VERSION=$(shell nvcc --version | tail -n 1 | cut -b 38-)
NVCC=nvcc -DCOMP='"NVCC $(NVCC_VERSION)"'

CUDALIBS=-L/usr/local/cuda/lib64/ -lcudart -lcusparse -lcurand

# List module object filenames
CUDA_OBJECTS= \
	     obj/cuda/config.o \
	     obj/cuda/data.o \
		  obj/cuda/dipole.o \
		  obj/cuda/exchange_fields.o \
		  obj/cuda/external_fields.o \
		  obj/cuda/finalize.o \
	     obj/cuda/initialize.o \
		  obj/cuda/internal.o \
		  obj/cuda/llg_heun.o \
	     obj/cuda/spin_fields.o \
	     obj/cuda/statistics.o

# The device linked object for CUDA
CUDA_DLINK= obj/cuda/vampire_cuda.o

# set precision for cuda calculation
CUDA_PR= -DCUDA_DP
#CUDA_PR= -DCUDA_SP
#CUDA_PR= -DCUDA_MP

CUDA_CFLAGS= -O3 -I./hdr/cusplibrary-0.5.1 -I./hdr -DCUDA -arch=sm_20 --use_fast_math -ftz=true

$(CUDA_OBJECTS): obj/%.o: src/%.cu
	$(NVCC) $< -c -o $@ $(CUDA_CFLAGS) $(CUDA_PR)

# Target for device link -- not working atm
#$(CUDA_DLINK): $(cuda_objects)
#	nvcc -arch=sm_30 -dlink $(cuda_objects) -o $@

# Target to link the CUDA enabled vampire
ifeq ($(shell uname -s), Darwin)
gcc-cuda: GCC_CFLAGS=$(LLVM_CFLAGS)
gcc-cuda: GCC_LDFLAGS=$(LLVM_LDFLAGS)
gcc-cuda: CUDALIBS=-lcudart -lcusparse -lcurand
endif
gcc-cuda: GCC_CFLAGS += -DCUDA
gcc-cuda: $(OBJECTS) $(CUDA_OBJECTS)
	$(NVCC) $(GCC_LDFLAGS) $(LIBS) ${CUDALIBS} $(OBJECTS) $(CUDA_OBJECTS) -o $(EXECUTABLE)

# CUDA debug target
gcc-cuda-debug: CUDA_CFLAGS += -DCUDA_DEBUG -DCUDA_SPIN_DEBUG
gcc-cuda-debug: gcc-cuda

# CUDA debug target
gcc-cuda-profile: CUDA_CFLAGS += -DCUDA_DEBUG -DCUDA_SPIN_DEBUG --profile
gcc-cuda-profile: gcc-cuda
