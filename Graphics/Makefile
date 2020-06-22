# module load gcc; module load plplot
# setenv SILO_HOME /usr/pppl/silo/4.10.2
#

F90 = gfortran
FOPTS = -ffree-line-length-none -O

PLPLOTLIB = -L$(PLPLOT_HOME)/lib -lplplotfortran

OBJ = bspline90_22_ver0.o init_value_run_plotter.o

silo.inc:
	sed /^C/d $(SILO_HOME)/include/silo.inc > silo.inc

bspline90_22_ver0.o: bspline90_22_ver0.f90
	$(F90) -c $(FOPTS) $<

init_value_run_plotter.o: init_value_run_plotter.f90 bspline90_22_ver0.o silo.inc
	$(F90) -c $(FOPTS) -I${PLPLOT_HOME}/lib/fortran/modules/plplot $<

xplot_init_value_run: $(OBJ)
	$(F90) $(FOPTS) -o $@ $(OBJ) $(PLPLOTLIB) -L$(SILO_HOME)/lib -lsilo

all: xplot_init_value_run

clean:
	rm -f *.o *.mod silo.inc

#ifort -c bspline90_22_ver0.f90
#ifort -c init_value_run_plotter.f90 -I. `PKG_CONFIG_PATH=/usr/local/plplot-5.11.1/lib/pkgconfig pkg-config --cflags --libs plplot-f95`
#ifort -O -o xplot_init_value_run init_value_run_plotter.o bspline90_22_ver0.o -I. `PKG_CONFIG_PATH=/usr/local/plplot-5.11.1/lib/pkgconfig pkg-config --cflags --libs plplot-f95` \
#  -L/Users/dsp/fortran_code_development/Silo-4.10.2/lib -lsilo
