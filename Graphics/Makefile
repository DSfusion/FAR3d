# PPPL CentOS 7 setup:
# module load gcc; module load plplot silo netcdf-c szip
# module load netcdf-fortran
#
# setenv LIBSTELLDIR /p/m3d/jbreslau/Stellarator/Github/STELLOPT/LIBSTELL/Release

F90 = gfortran
FOPTS = -ffree-line-length-none -O

PLPLOTLIB = -L$(PLPLOT_HOME)/lib -lplplotfortran

OBJ1 = bspline90_22_ver0.o init_value_run_plotter.o
OBJ2 = bspline90_22_ver0.o egn_plotter_jdqz_far3d.o
EXE = xplot_init_value_run xplot_egn_jdqz_far3d xbw_read

silo.inc:
	sed /^C/d $(SILO_DIR)/include/silo.inc > silo.inc

bspline90_22_ver0.o: bspline90_22_ver0.f90
	$(F90) -c $(FOPTS) $<

init_value_run_plotter.o: init_value_run_plotter.f90 bspline90_22_ver0.o silo.inc
	$(F90) -c $(FOPTS) -I${PLPLOT_HOME}/lib/fortran/modules/plplot $<

egn_plotter_jdqz_far3d.o: egn_plotter_jdqz_far3d.f90
	$(F90) -c $(FOPTS) -I${PLPLOT_HOME}/lib/fortran/modules/plplot $<

xplot_init_value_run: $(OBJ1)
	$(F90) $(FOPTS) -o $@ $(OBJ1) $(PLPLOTLIB) -L$(SILO_DIR)/lib -lsilo

xplot_egn_jdqz_far3d: $(OBJ2)
	$(F90) $(FOPTS) -o $@ $(OBJ2) $(PLPLOTLIB) -L$(SILO_DIR)/lib -lsilo

xbw_read: read_boozer_wout.f90
	$(F90) $(FOPTS) -o $@ $< -I$(LIBSTELLDIR) -L$(LIBSTELLDIR) -lstell \
	-L$(NETCDF_FORTRAN_HOME)/lib -lnetcdff

all: $(EXE)

clean:
	rm -f *.o *.mod silo.inc $(EXE)

#ifort -c bspline90_22_ver0.f90
#ifort -c init_value_run_plotter.f90 -I. `PKG_CONFIG_PATH=/usr/local/plplot-5.11.1/lib/pkgconfig pkg-config --cflags --libs plplot-f95`
#ifort -O -o xplot_init_value_run init_value_run_plotter.o bspline90_22_ver0.o -I. `PKG_CONFIG_PATH=/usr/local/plplot-5.11.1/lib/pkgconfig pkg-config --cflags --libs plplot-f95` \
#  -L/Users/dsp/fortran_code_development/Silo-4.10.2/lib -lsilo
