#-----------------------------------------------------------------------
# Makefile for FAR3d executable
# On PPPL cluster: module load intel
#-----------------------------------------------------------------------

FAR3D_HOME=..

#These are generic flags and commands.
#Comp = gfortran
#Flag = -ffree-line-length-none
Comp = ifort 
Flag = -free
Exc = xfar3d xfar3d_n
Opt = -O2

ifdef IMAS_PREFIX
 IMAS_FLAG = -DIMAS

 IMAS = -I$(IMAS_PREFIX)/include/$(Comp)
 I2MEX = -I$(I2MEX_HOME)/mod
 PSPLINE = -I$(PSPLINE_HOME)/include

 NTCC = -L$(NTCC_HOME)/lib -lold_xplasma -lxplasma2 -lfluxav -lr8bloat \
        -lsmlib -lezcdf -lcomput -lportlib
 LIBS = -lfftw3 -li2mex \
         $(NTCC) -L$(PSPLINE_HOME)/lib -lpspline -lnetcdff \
         -L$(IMAS_PREFIX)/lib -limas-$(Comp)

endif

#Library target
ARC = $(FAR3D_HOME)/lib/libfar3d.a

#Main program files
OBJS = Modules.o Main.o functions.o

#Setup subroutines files
OBJS2 = dfault.o setup.o setmod.o seteq.o etachi.o grid.o pert.o vmec.o ae_profiles.o inputlist.o inputlist_namelist.o

#Model subroutines
OBJS3 = b2lx.o block.o linstart.o lincheck.o linstep.o

#Operator subroutines
OBJS4 = clgam.o derivatives.o eqdy_dyeq.o om.o dlstar.o eigensolver_tools.o Laundamp_tools.o dlsq.o dlsq_r.o

#Output subroutines
OBJS5 = endrun.o energy.o output.o wrdump.o rddump.o

#Tool subroutines
OBJS6 = mult.o cnvt.o decbt.o solbt.o quadq.o mmlims.o numinc.o eqsplns.o elapsed_time.o fitter.o

libs: chkdirs $(ARC)

chkdirs:
	@test -d $(FAR3D_HOME)/lib || mkdir -p $(FAR3D_HOME)/lib

ifdef IMAS_PREFIX

#Transp2Far3D subroutines
OBJS7 = scxform.o far3d_wout_mod.o readtranspeq.o transp2i2mex.o ids_eqGet.o i2mex2far3d.o

scxform.o: $(TRANSP2FAR3D_HOME)/scxform.f90
	$(Comp) -assume nounderscore -cpp -c $(Opt) $(Flag) $<

far3d_wout_mod.o: $(TRANSP2FAR3D_HOME)/far3d_wout_mod.f90
	$(Comp) -cpp -c $(Opt) $(Flag) $<

readtranspeq.o: $(TRANSP2FAR3D_HOME)/readtranspeq.f90
	$(Comp) -cpp -c $(Opt) $(Flag) $<

transp2i2mex.o: $(TRANSP2FAR3D_HOME)/transp2i2mex.f90
	$(Comp) $(I2MEX) $(PSPLINE) -cpp -c $(Opt) $(Flag) $<

i2mex2far3d.o: $(TRANSP2FAR3D_HOME)/i2mex2far3d.f90 far3d_wout_mod.o
	$(Comp) $(I2MEX) -cpp -c $(Opt) $(Flag) $<

ids_eqGet.o: $(TRANSP2FAR3D_HOME)/ids_eqGet.f90
	$(Comp) $(IMAS) -cpp -c $(Opt) $(Flag) $<

endif

%.o: %.f90
	$(Comp) $(IMAS_FLAG) -cpp -c $(Opt) $(Flag) $<

$(ARC): $(OBJS) $(OBJS7) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6)
	ar -r $(ARC) $(OBJS) $(OBJS7) $(OBJS2) $(OBJS3) $(OBJS4) $(OBJS5) $(OBJS6)

#Code compilation
xfar3d: driver.f90 libs
	$(Comp) -o $@ $(Opt) $(Flag) -cpp $< -L$(FAR3D_HOME)/lib -lfar3d $(LIBS)

xfar3d_n: driver.f90 libs
	$(Comp) -o $@ $(Opt) $(Flag) -cpp -DNAMELIST_INPUT $< -L$(FAR3D_HOME)/lib -lfar3d $(LIBS)

all: $(Exc)

clean:
	rm -f *.o *.mod $(Exc) $(ARC)
