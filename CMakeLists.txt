cmake_minimum_required(VERSION 3.24)

option(FAR3D_ENABLE_GPU "Build the NVIDIA GPU implementation" OFF)
set(FAR3D_NVHPC_ROOT "" CACHE PATH
    "NVHPC installation root; normally discovered from PATH, NVHPC, or /opt/nvidia")
set(FAR3D_CUDA_VERSION "" CACHE STRING
    "CUDA Toolkit compatibility version passed to NVHPC (for example 13.1)")
set(FAR3D_GPU_ARCH "" CACHE STRING
    "NVHPC GPU architecture (for example cc80); empty uses the compiler default")

# Compiler selection has to happen before project() enables the languages.
# A GPU compatibility build with no explicit toolchain searches a standard
# NVHPC installation and uses its native compilers. MPI is attached below via
# the matching HPC-X wrapper, so the Makefile does not need compiler paths.
if(FAR3D_ENABLE_GPU AND NOT CMAKE_Fortran_COMPILER)
    set(_far3d_nvhpc_candidates)

    if(FAR3D_NVHPC_ROOT)
        list(APPEND _far3d_nvhpc_candidates "${FAR3D_NVHPC_ROOT}")
    else()
        if(DEFINED ENV{NVHPC} AND NOT "$ENV{NVHPC}" STREQUAL "")
            list(APPEND _far3d_nvhpc_candidates "$ENV{NVHPC}")
        endif()

        find_program(_far3d_nvfortran_in_path NAMES nvfortran NO_CACHE)
        if(_far3d_nvfortran_in_path)
            get_filename_component(_far3d_nvhpc_bin
                "${_far3d_nvfortran_in_path}" DIRECTORY)
            get_filename_component(_far3d_nvhpc_compilers
                "${_far3d_nvhpc_bin}" DIRECTORY)
            get_filename_component(_far3d_nvhpc_from_path
                "${_far3d_nvhpc_compilers}" DIRECTORY)
            list(APPEND _far3d_nvhpc_candidates "${_far3d_nvhpc_from_path}")
        endif()

        file(GLOB _far3d_nvhpc_in_opt LIST_DIRECTORIES TRUE
            "/opt/nvidia/hpc_sdk/*/*")
        list(SORT _far3d_nvhpc_in_opt COMPARE NATURAL ORDER DESCENDING)
        list(APPEND _far3d_nvhpc_candidates ${_far3d_nvhpc_in_opt})
    endif()
    list(REMOVE_DUPLICATES _far3d_nvhpc_candidates)

    foreach(_far3d_nvhpc_candidate IN LISTS _far3d_nvhpc_candidates)
        file(REAL_PATH "${_far3d_nvhpc_candidate}"
            _far3d_nvhpc_candidate_real)
        if(EXISTS "${_far3d_nvhpc_candidate}/compilers/bin/nvc" AND
           EXISTS "${_far3d_nvhpc_candidate}/compilers/bin/nvfortran")
            set(FAR3D_NVHPC_ROOT "${_far3d_nvhpc_candidate_real}" CACHE PATH
                "NVHPC installation root" FORCE)
            set(CMAKE_C_COMPILER
                "${FAR3D_NVHPC_ROOT}/compilers/bin/nvc" CACHE FILEPATH
                "C compiler" FORCE)
            set(CMAKE_Fortran_COMPILER
                "${FAR3D_NVHPC_ROOT}/compilers/bin/nvfortran" CACHE FILEPATH
                "Fortran compiler" FORCE)
            if(NOT CUDAToolkit_ROOT)
                if(FAR3D_CUDA_VERSION)
                    set(_far3d_cuda_candidates
                        "${FAR3D_NVHPC_ROOT}/cuda/${FAR3D_CUDA_VERSION}")
                else()
                    file(GLOB _far3d_cuda_candidates LIST_DIRECTORIES TRUE
                        "${FAR3D_NVHPC_ROOT}/cuda/*")
                    list(SORT _far3d_cuda_candidates COMPARE NATURAL
                        ORDER DESCENDING)
                endif()
                foreach(_far3d_cuda_candidate IN LISTS _far3d_cuda_candidates)
                    if(EXISTS "${_far3d_cuda_candidate}/bin/nvcc" AND
                       EXISTS "${_far3d_cuda_candidate}/include/cuda.h")
                        set(CUDAToolkit_ROOT "${_far3d_cuda_candidate}"
                            CACHE PATH "CUDA Toolkit root" FORCE)
                        break()
                    endif()
                endforeach()
            endif()
            break()
        endif()
    endforeach()
endif()

project(FAR3d VERSION 2.0 LANGUAGES C Fortran)

# Some CMake FindMPI versions do not include the modern Open MPI `mpifort`
# wrapper name in their candidate list. Locate a wrapper as a hint; FindMPI
# still interrogates it and validates the resulting imported target.
if(NOT MPI_Fortran_COMPILER)
    if(FAR3D_ENABLE_GPU AND FAR3D_NVHPC_ROOT AND
       EXISTS "${FAR3D_NVHPC_ROOT}/comm_libs/hpcx/bin/mpifort")
        set(MPI_Fortran_COMPILER
            "${FAR3D_NVHPC_ROOT}/comm_libs/hpcx/bin/mpifort" CACHE FILEPATH
            "MPI compiler wrapper used to discover the Fortran MPI interface")
    else()
        find_program(MPI_Fortran_COMPILER NAMES mpifort mpif90 mpif95
            DOC "MPI compiler wrapper used to discover the Fortran MPI interface")
    endif()
endif()

find_package(MPI REQUIRED COMPONENTS Fortran)
find_package(OpenMP REQUIRED COMPONENTS Fortran)
find_package(OpenACC REQUIRED COMPONENTS Fortran)

set(far3d_sources
    src/MPI_tools/transfer.f90
    src/bigf.f90
    src/block_mod.f90
    src/cnvt.f90
    src/dbyd.f90
    src/dbydb.f90
    src/delstar.f90
    src/dfault.f90
    src/diagnostics.f90
    src/divv.f90
    src/dlsq.f90
    src/dlstar.f90
    src/equilibrium.f90
    src/far3d.f90
    src/findr.f90
    src/findrf.f90
    src/fit_mod.f90
    src/gam_mod.f90
    src/globals.f90
    src/grid.f90
    src/initialize.f90
    src/lincheck.f90
    src/linstart.f90
    src/mult_mod.f90
    src/multb_mod.f90
    src/om_mod.f90
    src/output_mod.f90
    src/scratch.f90
    src/solve.f90
    src/tools.f90)

if(FAR3D_ENABLE_GPU)
    if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(NVHPC|PGI)$")
        message(FATAL_ERROR
            "FAR3D_ENABLE_GPU requires the NVHPC Fortran compiler; "
            "the selected compiler is ${CMAKE_Fortran_COMPILER_ID} "
            "(${CMAKE_Fortran_COMPILER}).")
    endif()

    if(FAR3D_CUDA_VERSION)
        find_package(CUDAToolkit ${FAR3D_CUDA_VERSION} EXACT REQUIRED)
    else()
        find_package(CUDAToolkit REQUIRED)
    endif()
    list(APPEND far3d_sources src/matrix_gpu.f90)
else()
    list(APPEND far3d_sources src/matrix_cpu.f90)
endif()

add_executable(far3d ${far3d_sources})
set_target_properties(far3d PROPERTIES
    OUTPUT_NAME far3d.x
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules")

target_include_directories(far3d PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/modules")
target_link_libraries(far3d PRIVATE
    MPI::MPI_Fortran
    OpenMP::OpenMP_Fortran
    OpenACC::OpenACC_Fortran)

# Runtime diagnostics for the debug configuration. Optimization and debug
# symbols themselves are supplied by CMake's standard build-type flags.
target_compile_options(far3d PRIVATE
    "$<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:GNU>>:-fcheck=all;-fbacktrace>"
    "$<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:NVHPC,PGI>>:-Mbounds;-Mchkptr;-Mchkstk;-traceback>")

if(FAR3D_ENABLE_GPU)
    target_link_libraries(far3d PRIVATE CUDA::cudart)

    set(far3d_gpu_options -acc=gpu -gpu=mem:separate)
    if(FAR3D_CUDA_VERSION)
        list(APPEND far3d_gpu_options -gpu=cuda${FAR3D_CUDA_VERSION})
    endif()
    if(FAR3D_GPU_ARCH)
        list(APPEND far3d_gpu_options -gpu=${FAR3D_GPU_ARCH})
    endif()

    # -cuda enables CUDA Fortran (the GPU matrix uses the cudafor module) and
    # must be present at both compilation and final link time.
    target_compile_options(far3d PRIVATE ${far3d_gpu_options} -cuda)
    target_link_options(far3d PRIVATE ${far3d_gpu_options} -cuda)
endif()

install(TARGETS far3d RUNTIME DESTINATION bin)

if(FAR3D_ENABLE_GPU)
    set(far3d_variant GPU)
else()
    set(far3d_variant CPU)
endif()
message(STATUS "FAR3D variant: ${far3d_variant}")
message(STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER_ID} ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "MPI Fortran compiler: ${MPI_Fortran_COMPILER}")
