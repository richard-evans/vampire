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
	obj/cuda/hamr_fields.o \
	obj/cuda/initialize.o \
	obj/cuda/internal.o \
	obj/cuda/llg_heun.o \
	obj/cuda/monte_carlo.o \
	obj/cuda/spin_fields.o \
	obj/cuda/statistics.o \
	obj/cuda/thermal_fields.o \
	obj/cuda/transfer.o 

CUDA_ARCH=\
-gencode arch=compute_35,code=sm_35 \
-gencode arch=compute_50,code=sm_50 \
-gencode arch=compute_50,code=compute_50 \
-gencode arch=compute_61,code=sm_61 \
-gencode arch=compute_61,code=compute_61 \
-gencode arch=compute_70,code=compute_70 \
-gencode arch=compute_75,code=compute_75

# The device linked object for CUDA
CUDA_DLINK= obj/cuda/vampire_cuda.o

# set precision for cuda calculation
CUDA_PR= -DCUDA_DP
#CUDA_PR= -DCUDA_SP
#CUDA_PR= -DCUDA_MP

# set matrix format
CUDA_MATRIX_TYPE = CSR
#CUDA_MATRIX_TYPE = DIA
#CUDA_MATRIX_TYPE = ELL

#CUDA_CFLAGS= -O3 -I./hdr/cusplibrary-0.5.1 -I./hdr -DCUDA $(CUDA_ARCH) --use_fast_math -ftz=true -std=c++11

#CUDA 8.0+
#CUDA_CFLAGS= -O3 -I./hdr -DCUDA $(CUDA_ARCH) --use_fast_math -ftz=true -std=c++11 -Wno-deprecated-gpu-targets
#CUDA 10.0 +
CUDA_CFLAGS= -O3 -I./hdr -DCUDA $(CUDA_ARCH) --use_fast_math -ftz=true -std=c++14 -Wno-deprecated-gpu-targets

$(CUDA_OBJECTS): obj/%.o: src/%.cu
	$(NVCC) $< -dc -o $@ $(CUDA_CFLAGS) $(CUDA_PR) -DCUDA_MATRIX=$(CUDA_MATRIX_TYPE)

# Target to link the CUDA enabled vampire
ifeq ($(shell uname -s), Darwin)
cuda: GCC_CFLAGS=$(LLVM_CFLAGS)
cuda: GCC_LDFLAGS=$(LLVM_LDFLAGS)
cuda: CUDALIBS=-lcudart -lcusparse -lcurand
endif
cuda: GCC_CFLAGS += -DCUDA
cuda: $(OBJECTS) $(CUDA_OBJECTS)
	$(NVCC) $(CUDA_CFLAGS) $(GCC_LDFLAGS) $(LIBS) ${CUDALIBS} $(OBJECTS) $(CUDA_OBJECTS) -o $(CUDAEXECUTABLE)

# CUDA debug target
cuda-debug: CUDA_CFLAGS += -DCUDA_DEBUG -DCUDA_SPIN_DEBUG
cuda-debug: cuda

# CUDA debug target
cuda-profile: CUDA_CFLAGS += -DCUDA_DEBUG -DCUDA_SPIN_DEBUG --profile
cuda-profile: cuda
