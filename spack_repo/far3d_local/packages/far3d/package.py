# Copyright (c) 2024 FAR3d developers
#
# SPDX-License-Identifier: MIT

from spack_repo.builtin.build_systems.cmake import CMakePackage

from spack.package import depends_on, license, requires, variant, version


class Far3d(CMakePackage):
    """Parallel gyrofluid code for nonlinear simulations of energetic
    particle-driven instabilities in three-dimensional configurations."""

    homepage = "https://github.com/DSfusion/FAR3d"
    git = "https://github.com/DSfusion/FAR3d.git"

    license("MIT")

    # This version is normally associated with the already-cloned checkout by
    # `spack develop --no-clone`; the Git metadata is only a fallback fetcher.
    version("develop", branch="master")

    variant("cuda", default=False, description="Build the NVHPC CUDA Fortran implementation")

    depends_on("c", type="build")
    depends_on("fortran", type="build")
    depends_on("cmake@3.24:", type="build")
    depends_on("mpi")
    depends_on("cuda", when="+cuda", type=("build", "link", "run"))

    requires("%nvhpc", when="+cuda", msg="far3d+cuda requires CUDA Fortran from NVHPC")

    def cmake_args(self):
        spec = self.spec
        args = [
            # Use the MPI wrappers as the project compilers. Spack guarantees
            # that this MPI dependency is paired with the selected compiler.
            self.define("CMAKE_C_COMPILER", spec["mpi"].mpicc),
            self.define("CMAKE_Fortran_COMPILER", spec["mpi"].mpifc),
            self.define_from_variant("FAR3D_ENABLE_GPU", "cuda"),
        ]

        if "+cuda" in spec:
            cuda_version = str(spec["cuda"].version.up_to(2))
            args.extend(
                [
                    self.define("CUDAToolkit_ROOT", spec["cuda"].prefix),
                    self.define("FAR3D_CUDA_VERSION", cuda_version),
                ]
            )

        return args
