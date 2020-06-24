# PPPL cluster setup:
#  module load intel

FC     = ifort
FFLAGS = -fixed
FOPTS  = -O2
MKL    = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core \
         -lmkl_sequential -lpthread -lm

MAIN   = grandom_mod.o symtrd_mod.o myjdqz_mod_TAEFL_cmplx.o \
         main_jdqz_TAEFL_cmplx.o

EXE    = xEigen

%.o: %.f90
	$(FC) $(FFLAGS) $(FOPTS) -c $<

xEigen: $(MAIN)
	cd LIB_JDQZ; make libjdqz.a
	$(FC) $(FOPTS) -o $@ $(MAIN) -LLIB_JDQZ -ljdqz $(MKL)

all: $(EXE)

clean:
	cd LIB_JDQZ; make clean
	rm -f *.o *.mod $(EXE)
