# PPPL cluster setup:
#  module load intel

FAR3D_HOME=..

FC     = ifort
FFLAGS = -fixed
FOPTS  = -O2
MKL    = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core \
         -lmkl_sequential -lpthread -lm

ARC    = $(FAR3D_HOME)/lib/libfar3deigen.a

OBJS   = grandom_mod.o symtrd_mod.o myjdqz_mod_TAEFL_cmplx.o \
         main_jdqz_TAEFL_cmplx.o

EXE    = xEigen

all: $(EXE)

libs: chkdirs $(ARC)
	cd LIB_JDQZ; make libjdqz.a

chkdirs:
	@test -d $(FAR3D_HOME)/lib || mkdir -p $(FAR3D_HOME)/lib

$(ARC): $(OBJS)
	ar -r $(ARC) $(OBJS)

%.o: %.f90
	$(FC) $(FFLAGS) $(FOPTS) -c $<

xEigen: driver.f90 libs
	$(FC) $(FOPTS) -o $@ $< -L$(FAR3D_HOME)/lib -lfar3deigen \
	 -LLIB_JDQZ -ljdqz $(MKL)

refresh: clean all

clean:
	cd LIB_JDQZ; make clean
	rm -f *.o *.mod $(EXE) $(ARC)
