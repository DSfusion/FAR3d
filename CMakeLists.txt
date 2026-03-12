cmake_minimum_required(VERSION 3.10)

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
# Source files
# -------------------------
# Automatically include all Fortran source files
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.f90 src/*.F90 src/*.f)

# -------------------------
# Add executable
# -------------------------
add_executable(far3d.x ${SOURCES})
target_include_directories(far3d.x PRIVATE)
# -------------------------
# Link libraries
# -------------------------
target_link_libraries(far3d.x
    MPI::MPI_Fortran
    OpenMP::OpenMP_Fortran
)
# -------------------------
# Optional: Compiler flags
# -------------------------
# You can use CMake flags set from the outside (e.g., via Makefile), or set defaults here.
# target_compile_options(far3d.x PRIVATE -g -O3 -fopenmp)

message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_Fortran_FLAGS = ${CMAKE_Fortran_FLAGS}")
message(STATUS "CMAKE_Fortran_FLAGS_RELEASE = ${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "CMAKE_Fortran_FLAGS_DEBUG = ${CMAKE_Fortran_FLAGS_DEBUG}")
