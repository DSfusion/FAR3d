# On PPPL cluster (CentOS 6 Linux, csh):
# module load intel
# setenv STELLOPT_PATH /p/stellopt/SOURCES/stellopt250_intel
# module load lapack netcdf-c
# module load netcdf-fortran

FC    = ifort
FOPT = -O2 -fixed
LIBSTELLDIR = $(STELLOPT_PATH)/LIBSTELL/Release

MOD = booz_params.o booz_persistent.o
OBJ = read_wout_booz.o boozer_metric.o write_boozmn.o write_boozmn_cdf.o free_mem_boozer.o \
      boozer_coords.o allocate_boozer.o write_polcut.o trigfunc.o vcoords_gijb.o vcoords.o \
      boozer_gij.o transpmn.o harfun.o booz_jac.o boozer.o setup_booz.o foranl.o root.o
LIBS =  -L. -lstell -L${NETCDF_FORTRAN_HOME}/lib -lnetcdff -L${LAPACK_HOME}/lib -llapack -lrefblas

EXE = xbooz_xform

xbooz_xform: boozer_xform.f $(MOD) $(OBJ) libstell.a
	$(FC) $(FOPT) -I. -I$(LIBSTELLDIR) -o $@ $< $(MOD) $(OBJ) $(LIBS)

%.o: %.f
	$(FC) $(FOPT) -I$(LIBSTELLDIR) -c $<

libstell.a:
	$(FC) -O -c $(STELLOPT_PATH)/LIBSTELL/Sources/Ezcdf/handle_err.f90
	ar -r libstell.a $(LIBSTELLDIR)/*.o handle_err.o

all: $(EXE)

clean:
	rm -f *.o *.mod libstell.a $(EXE)
