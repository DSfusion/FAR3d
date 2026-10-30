MAKE_CMD="make"
ifdef VERBOSE
        MAKE_CMD+=" VERBOSE=1"
endif

# BLAS/LAPACK backend, e.g. make LAPACK=None, LAPACK=Accelerate, LAPACK=MKL, LAPACK=Reference or LAPACK=NPC
CMAKE_OPTS=
ifdef LAPACK
        CMAKE_OPTS+=-DFAR3D_LAPACK=$(LAPACK)
endif

all: ../build
	@cd ../build; \
		FC=$(F90) cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_Fortran_COMPILER=mpif90 $(CMAKE_OPTS) ..; \
		$(MAKE_CMD)

debug: ../build
	@cd ../build; \
		FC=$(F90) cmake -DCMAKE_BUILD_TYPE=Debug $(CMAKE_OPTS) ..; \
		$(MAKE_CMD)

../build:
	@mkdir ../build;

clean:  
	@echo "Cleaning build folder"
	@rm -rf ../build
	@rm -f far3d.x
