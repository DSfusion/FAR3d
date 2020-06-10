#-----------------------------------------------------------------------
# Makefile for FAR3d executable
# On PPPL cluster: module load intel
#-----------------------------------------------------------------------

#These are generic flags and commands.
#Comp="gfortran" 
#Flag=" -ffree-line-length-none"
Comp = ifort 
Flag = -free
Exc = xfar3d 
Opt = -O2

#Main program files
OBJS = Modules.o Main.o

#Setup subroutines files
OBJS2 = inputlist.o dfault.o setup.o setmod.o seteq.o etachi.o grid.o pert.o vmec.o ae_profiles.o

#Model subroutines
OBJS3 = b2lx.o block.o linstart.o lincheck.o linstep.o

#Operator subroutines
OBJS4 = clgam.o derivatives.o eqdy_dyeq.o om.o dlstar.o eigensolver_tools.o Laundamp_tools.o dlsq.o

#Output subroutines
OBJS5 = endrun.o energy.o output.o wrdump.o rddump.o

#Tool subroutines
OBJS6 = mult.o cnvt.o decbt.o solbt.o quadq.o functions.o mmlims.o numinc.o eqsplns.o elapsed_time.o fitter.o

%.o: %.f90
	$(Comp) -c $(Opt) $(Flag) $<

#Code compilation
$(Exc): $(OBJS) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6)
	$(Comp) -o $@ $(Opt) $(Flag) $(OBJS) $(OBJS2) $(OBJS3) $(OBJS4) \
	$(OBJS5) $(OBJS6)

clean:
	rm -f *.o *.mod xfar3d
