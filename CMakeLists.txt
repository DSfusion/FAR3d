cmake_minimum_required(VERSION 3.13)

project(far3d LANGUAGES Fortran C)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(DEBUG_FLAGS "-g -fcheck=all -fbacktrace")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(DEBUG_FLAGS "-g -check all -traceback")
else()
    message(WARNING "Unknown Fortran compiler, using default debug flags")
    set(DEBUG_FLAGS "-g")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${DEBUG_FLAGS}")

# -------------------------
# Find MPI
# -------------------------
find_package(MPI REQUIRED)

# -------------------------
# Find OpenMP (optional, but often used with MPI)
# -------------------------
find_package(OpenMP REQUIRED)

# -------------------------
# BLAS / LAPACK
# -------------------------
# FAR3D_LAPACK selects the BLAS/LAPACK backend for the factorisation (decbt) in
# matrix.f90; the solves in the time step (solbt_all) do not use BLAS:
#   None       - built-in dec/sol routines, no BLAS/LAPACK
#   Accelerate - Apple Accelerate framework (default on macOS)
#   MKL        - Intel oneMKL (default elsewhere), found via its CMake package
#                (MKL_DIR, MKLROOT or MKL_ROOT, e.g. after sourcing oneAPI setvars.sh)
#   Reference  - reference BLAS/LAPACK from FAR3D_REFLAPACK_ROOT
#   NPC        - LAPACK from NAG Portable Cores (libnpc_lp64 in FAR3D_NPC_ROOT/lib),
#                BLAS from sequential MKL (the BLAS backend NPC was built with)
if(APPLE)
    set(FAR3D_LAPACK_DEFAULT "Accelerate")
else()
    set(FAR3D_LAPACK_DEFAULT "MKL")
endif()
set(FAR3D_LAPACK "${FAR3D_LAPACK_DEFAULT}" CACHE STRING "BLAS/LAPACK backend: None, Accelerate, MKL, Reference or NPC")
set_property(CACHE FAR3D_LAPACK PROPERTY STRINGS None Accelerate MKL Reference NPC)
set(FAR3D_REFLAPACK_ROOT "/opt/homebrew/opt/lapack" CACHE PATH "Install prefix of the reference BLAS/LAPACK")
set(FAR3D_NPC_ROOT "" CACHE PATH "NPC build or install directory containing lib/libnpc_lp64")
# libnpc_lp64.so depends on the LLVM OpenMP runtime (libomp)
set(FAR3D_NPC_OMP_DIR "" CACHE PATH "Directory containing libomp.so needed by NPC")
# gnu_thread shares libgomp with the OpenMP runtime used by gfortran
set(FAR3D_MKL_THREADING "gnu_thread" CACHE STRING "MKL threading layer: gnu_thread, intel_thread or sequential")
set_property(CACHE FAR3D_MKL_THREADING PROPERTY STRINGS gnu_thread intel_thread sequential)

set(LAPACK_LIBS "")
if(FAR3D_LAPACK STREQUAL "None")
elseif(FAR3D_LAPACK STREQUAL "Accelerate")
    find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
    set(LAPACK_LIBS ${ACCELERATE_FRAMEWORK})
elseif(FAR3D_LAPACK STREQUAL "Reference")
    find_library(REF_LAPACK_LIB lapack PATHS ${FAR3D_REFLAPACK_ROOT}/lib NO_DEFAULT_PATH REQUIRED)
    find_library(REF_BLAS_LIB blas PATHS ${FAR3D_REFLAPACK_ROOT}/lib NO_DEFAULT_PATH REQUIRED)
    set(LAPACK_LIBS ${REF_LAPACK_LIB} ${REF_BLAS_LIB})
elseif(FAR3D_LAPACK STREQUAL "MKL")
    # LP64 interface: matrix.f90 passes default (32-bit) integers to BLAS/LAPACK
    set(MKL_INTERFACE "lp64")
    set(MKL_THREADING "${FAR3D_MKL_THREADING}")
    find_package(MKL CONFIG REQUIRED HINTS $ENV{MKLROOT})
    set(LAPACK_LIBS MKL::MKL)
elseif(FAR3D_LAPACK STREQUAL "NPC")
    find_library(NPC_LIB npc_lp64 PATHS ${FAR3D_NPC_ROOT}/lib NO_DEFAULT_PATH REQUIRED)
    find_library(NPC_OMP_LIB omp PATHS ${FAR3D_NPC_OMP_DIR} REQUIRED)
    # dgemm/dgemv come from MKL; NPC itself links mkl_sequential, so match it
    set(MKL_INTERFACE "lp64")
    set(MKL_THREADING "sequential")
    find_package(MKL CONFIG REQUIRED HINTS $ENV{MKLROOT})
    # NPC first so that its dgetrf/dgetrs take precedence over MKL's
    set(LAPACK_LIBS ${NPC_LIB} MKL::MKL ${NPC_OMP_LIB})
else()
    message(FATAL_ERROR "Unknown FAR3D_LAPACK='${FAR3D_LAPACK}' (use None, Accelerate, MKL, Reference or NPC)")
endif()
message(STATUS "FAR3D_LAPACK = ${FAR3D_LAPACK} ${LAPACK_LIBS}")

# -------------------------
# Source files
# -------------------------
# Automatically include all Fortran source files
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.f90 src/*.F90 src/*.f)

# -------------------------
# Add executable
# -------------------------
add_executable(far3d.x ${SOURCES})
target_include_directories(far3d.x PRIVATE)
# matrix.f90 selects the LAPACK or built-in factorisation with the preprocessor
set_source_files_properties(src/matrix.f90 PROPERTIES Fortran_PREPROCESS ON)
if(NOT FAR3D_LAPACK STREQUAL "None")
    target_compile_definitions(far3d.x PRIVATE FAR3D_USE_LAPACK)
endif()
# -------------------------
# Link libraries
# -------------------------
target_link_libraries(far3d.x
    MPI::MPI_Fortran
    OpenMP::OpenMP_Fortran
    ${LAPACK_LIBS}
)
# -------------------------
# Optional: Compiler flags
# -------------------------
# You can use CMake flags set from the outside (e.g., via Makefile), or set defaults here.
# target_compile_options(far3d.x PRIVATE -g -O3)

message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE = ${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG = ${CMAKE_Fortran_FLAGS_DEBUG}")

# -------------------------
# CPU-specific code generation
# -------------------------
# Generate code for the CPU of the build machine (e.g. AVX-512); the executable
# may not run on older CPUs
option(FAR3D_MARCH_NATIVE "Optimise for the CPU of the build machine (-march=native / -xHost)" OFF)
if(FAR3D_MARCH_NATIVE)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        target_compile_options(far3d.x PRIVATE -xHost)
    else()
        target_compile_options(far3d.x PRIVATE -march=native)
    endif()
endif()
message(STATUS "FAR3D_MARCH_NATIVE = ${FAR3D_MARCH_NATIVE}")
