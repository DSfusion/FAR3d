FAR3d version 2.0

(This is a parallel version that requires MPI)

The earlier version of FAR3d has now been moved to the PreviousFar3d branch

Both CPU and GPU versions are included in this code.

## Build

FAR3d requires CMake 3.24 or newer, a Fortran compiler, MPI, and OpenMP.
The GPU variant additionally requires NVHPC, OpenACC, and CUDA. Configure from
the repository root; CMake discovers the MPI and parallel-runtime flags instead
of using an MPI compiler wrapper as the project compiler.

The historical entry points remain available from `src`:

```sh
make -C src all      # Release CPU build in build-release
make -C src debug    # Debug CPU build in build-debug
make -C src gpu      # Release NVIDIA GPU build in build-gpu
```

For the bare `gpu` target, CMake discovers NVHPC from `nvfortran` in `PATH`,
the `NVHPC` environment variable, or a versioned installation below
`/opt/nvidia/hpc_sdk`. For another location, pass
`CMAKE_CONFIGURE_ARGS="-DFAR3D_NVHPC_ROOT=/path/to/nvhpc"` to Make.

The CPU build selects `src/matrix_cpu.f90`. The GPU build selects
`src/matrix_gpu.f90`, requires NVHPC because it uses CUDA Fortran, and can be
configured directly with `-DFAR3D_ENABLE_GPU=ON`.

The separate `runconf.cpu` and `runconf.gpu` scripts configure Release builds
from their respective build directories. To configure and build the CPU
version with GCC/GFortran through the system MPI wrappers:

```sh
mkdir build-cpu
cd build-cpu
. ../runconf.cpu
make -j 8
```

To configure and build the GPU version with NVHPC 26.3 and CUDA 13.1:

```sh
mkdir build-gpu
cd build-gpu
. ../runconf.gpu
make -j 8
```

The resulting executable is `far3d.x` in the selected build directory. The
scripts contain the compiler, MPI, NVHPC, and CUDA paths for the installations
used by this checkout; edit those paths when using a different installation.

## Spack build from this checkout

The bundled `spack_repo` contains a local FAR3d package definition. Register
it once with Spack, then create an environment that associates the package's
`develop` version with this clone:

```sh
spack repo add "$PWD/spack_repo/far3d_local"
spack env create -d .spack-env
spack -e .spack-env develop --no-clone --path "$PWD" far3d@develop
```

Add and install either the CPU variant:

```sh
spack -e .spack-env add far3d@develop~cuda
spack -e .spack-env concretize
spack -e .spack-env install
```

or the CUDA Fortran variant:

```sh
spack -e .spack-env add far3d@develop+cuda
spack -e .spack-env concretize
spack -e .spack-env install
```

`far3d+cuda` requires an NVHPC compiler and a CUDA package known to Spack. For
example, select them and an NVIDIA compute capability with
`far3d@develop+cuda cuda_arch=80 %nvhpc ^cuda@13.1`. The recipe uses the
selected MPI package's `mpicc` and `mpifc` wrappers as the CMake compilers, so
MPI and the underlying Fortran compiler remain ABI-compatible.
