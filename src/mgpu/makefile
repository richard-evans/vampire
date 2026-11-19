#--------------------------------------------------------------
#      Makefile for the multi-GPU (MPI + CUDA) version
#--------------------------------------------------------------
# make cuda-multigpu builds vampire-cuda-multigpu. All sources are compiled
# again with MPI and CUDA enabled into their own objects (suffix _mgpu.o), so
# the serial, parallel and cuda versions are not affected. The MPI compiler
# wrapper (mpicxx) is used as host compiler for nvcc.

MGPUEXECUTABLE=vampire-cuda-multigpu
MGPU_NVCC=$(NVCC) -ccbin mpicxx -DMPICF -DCUDA_MULTIGPU

MGPU_HOST_OBJECTS=$(OBJECTS:.o=_mgpu.o) \
	obj/mgpu/data_mgpu.o \
	obj/mgpu/halo_mgpu.o \
	obj/mgpu/initialize_mgpu.o

MGPU_CUDA_OBJECTS=$(CUDA_OBJECTS:.o=_mgpu.o) \
	obj/mgpu/kernels_mgpu.o \
	obj/mgpu/thermal_mgpu.o

$(MGPU_HOST_OBJECTS): obj/%_mgpu.o: src/%.cpp
	$(MPICC) -c -o $@ $(GCC_CFLAGS) -DCUDA -DCUDA_MULTIGPU $(OPTIONS) $<

$(MGPU_CUDA_OBJECTS): obj/%_mgpu.o: src/%.cu
	$(MGPU_NVCC) $< -dc -o $@ $(CUDA_CFLAGS) $(CUDA_PR) -DCUDA_MATRIX=$(CUDA_MATRIX_TYPE)

cuda-multigpu: $(MGPU_HOST_OBJECTS) $(MGPU_CUDA_OBJECTS)
	$(MGPU_NVCC) $(CUDA_CFLAGS) $(GCC_LDFLAGS) $(LIBS) $(CUDALIBS) $(MGPU_HOST_OBJECTS) $(MGPU_CUDA_OBJECTS) -o $(MGPUEXECUTABLE)
