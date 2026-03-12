MAKE_CMD := make
ifdef VERBOSE
  MAKE_CMD += VERBOSE=1
endif

# Build directories
BUILD_RELEASE := ../build-release
BUILD_DEBUG   := ../build-debug
BUILD_GPU     := ../build-gpu

# GPU build flags
GPU_FC      := ftn
GPU_FFLAGS  := -O3 -mp -acc=gpu -gpu=managed -Minstrument -traceback -cuda
GPU_LDFLAGS := -lnvhpcwrapnvtx -cuda

.PHONY: all debug gpu clean

all: $(BUILD_RELEASE)
	@cd $(BUILD_RELEASE) && \
	cmake -DCMAKE_BUILD_TYPE=Release \
	  -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx \
	  -DCMAKE_Fortran_COMPILER=mpif90 .. && \
	$(MAKE_CMD)

debug: $(BUILD_DEBUG)
	@cd $(BUILD_DEBUG) && \
	cmake -DCMAKE_BUILD_TYPE=Debug .. && \
	$(MAKE_CMD)

gpu: $(BUILD_GPU)
	@cd $(BUILD_GPU) && \
	export OMPI_MCA_opal_cuda_support=0 && \
	cmake -DCMAKE_BUILD_TYPE=Release \
	  -DCMAKE_Fortran_COMPILER=$(GPU_FC) \
	  -DCMAKE_Fortran_FLAGS="$(GPU_FFLAGS)" \
	  -DCMAKE_EXE_LINKER_FLAGS="$(GPU_LDFLAGS)" \
	  .. && \
	$(MAKE_CMD)

# Directory creation rules
$(BUILD_RELEASE) $(BUILD_DEBUG) $(BUILD_GPU):
	@mkdir -p $@

clean:
	@echo "Cleaning build folders"
	@rm -rf $(BUILD_RELEASE) $(BUILD_DEBUG) $(BUILD_GPU)
	@rm -f far3d.x

