# Compatibility entry points for the historical "make all/debug/gpu" build.
# Compiler, MPI, OpenMP, OpenACC, and CUDA discovery is owned by CMake.

CMAKE ?= cmake
CMAKE_CONFIGURE_ARGS ?=
BUILD_PARALLEL_LEVEL ?=

ROOT_DIR := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)
BUILD_RELEASE := $(ROOT_DIR)/build-release
BUILD_DEBUG := $(ROOT_DIR)/build-debug
BUILD_GPU := $(ROOT_DIR)/build-gpu

ifneq ($(strip $(BUILD_PARALLEL_LEVEL)),)
  BUILD_PARALLEL_ARG := --parallel $(BUILD_PARALLEL_LEVEL)
else
  BUILD_PARALLEL_ARG := --parallel
endif

.PHONY: all debug gpu clean

all:
	$(CMAKE) $(CMAKE_CONFIGURE_ARGS) -S "$(ROOT_DIR)" -B "$(BUILD_RELEASE)" \
	  -DCMAKE_BUILD_TYPE=Release -DFAR3D_ENABLE_GPU=OFF
	$(CMAKE) --build "$(BUILD_RELEASE)" $(BUILD_PARALLEL_ARG)

debug:
	$(CMAKE) $(CMAKE_CONFIGURE_ARGS) -S "$(ROOT_DIR)" -B "$(BUILD_DEBUG)" \
	  -DCMAKE_BUILD_TYPE=Debug -DFAR3D_ENABLE_GPU=OFF
	$(CMAKE) --build "$(BUILD_DEBUG)" $(BUILD_PARALLEL_ARG)

gpu:
	$(CMAKE) --fresh $(CMAKE_CONFIGURE_ARGS) -S "$(ROOT_DIR)" -B "$(BUILD_GPU)" \
	  -DCMAKE_BUILD_TYPE=Release -DFAR3D_ENABLE_GPU=ON
	$(CMAKE) --build "$(BUILD_GPU)" $(BUILD_PARALLEL_ARG)

clean:
	$(CMAKE) -E rm -rf "$(BUILD_RELEASE)" "$(BUILD_DEBUG)" "$(BUILD_GPU)"
