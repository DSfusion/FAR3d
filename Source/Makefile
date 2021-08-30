#-----------------------------------------------------------------------
# Makefile for FAR3d executable
# On PPPL cluster: module load intel
#-----------------------------------------------------------------------

FAR3D_HOME=..

#These are generic flags and commands.
#Comp="gfortran" 
#Flag=" -ffree-line-length-none"
Comp = ifort 
Flag = -free
Exc = xfar3d xfar3d_n
Opt = -O2

#Library target
ARC = $(FAR3D_HOME)/lib/libfar3d.a

#Main program files
OBJS = Modules.o Main.o

#Setup subroutines files
OBJS2 = dfault.o setup.o setmod.o seteq.o etachi.o grid.o pert.o vmec.o ae_profiles.o inputlist.o inputlist_namelist.o

#Model subroutines
OBJS3 = b2lx.o block.o linstart.o lincheck.o linstep.o

#Operator subroutines
OBJS4 = clgam.o derivatives.o eqdy_dyeq.o om.o dlstar.o eigensolver_tools.o Laundamp_tools.o dlsq.o dlsq_r.o

#Output subroutines
OBJS5 = endrun.o energy.o output.o wrdump.o rddump.o

#Tool subroutines
OBJS6 = mult.o cnvt.o decbt.o solbt.o quadq.o functions.o mmlims.o numinc.o eqsplns.o elapsed_time.o fitter.o

libs: chkdirs $(ARC)

chkdirs:
	@test -d $(FAR3D_HOME)/lib || mkdir -p $(FAR3D_HOME)/lib

%.o: %.f90
	$(Comp) -c $(Opt) $(Flag) $<

$(ARC): $(OBJS) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6)
	ar -r $(ARC) $(OBJS) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6)

#Code compilation
xfar3d: driver.f90 libs
	$(Comp) -o $@ $(Opt) $(Flag) -cpp $< -L$(FAR3D_HOME)/lib -lfar3d

xfar3d_n: driver.f90 libs
	$(Comp) -o $@ $(Opt) $(Flag) -cpp -DNAMELIST_INPUT $< -L$(FAR3D_HOME)/lib -lfar3d

all: $(Exc)

clean:
	rm -f *.o *.mod $(Exc) $(ARC)
