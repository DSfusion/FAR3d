#!/bin/bash
#-----------------------------------------------------------------------
# Makefile for FAR3d executable
#-----------------------------------------------------------------------

#These are generic flags and commands.
Comp="gfortran" 
Exc=" xfar3d" 
Opt=" -O2"
Flag1=" -ffree-line-length-none -g -frange-check -fbounds-check"
Flag2=" -fbacktrace -ffpe-trap=invalid -ftrapv -ffree-line-length-none"

#Main program files
OBJS=" Modules.f90 Main.f90" 
#Setup subroutines files
OBJS2=" inputlist.f90 dfault.f90 setup.f90 setmod.f90 seteq.f90 etachi.f90 grid.f90 pert.f90 vmec.f90 ae_profiles.f90"
#Model subroutines
OBJS3=" b2lx.f90 block.f90 linstart.f90 lincheck.f90 linstep.f90"
#Operator subroutines
OBJS4=" clgam.f90 derivatives.f90 eqdy_dyeq.f90 om.f90 dlstar.f90 eigensolver_tools.f90 Laundamp_tools.f90 dlsq.f90 dlsq_r.f90"
#Output subroutines
OBJS5=" endrun.f90 energy.f90 output.f90 wrdump.f90 rddump.f90"
#Tool subroutines
OBJS6=" mult.f90 cnvt.f90 decbt.f90 solbt.f90 quadq.f90 functions.f90 mmlims.f90 numinc.f90 eqsplns.f90 elapsed_time.f90 fitter.f90"

#Code compilation
$Comp -o $Exc $Opt $Flag1 $Flag2 $OBJS $OBJS2 $OBJS3 $OBJS4 $OBJS5 $OBJS6       

#-----------------------------------------------------------------------
#                             Dependencies
#-----------------------------------------------------------------------
#Main.o:                        Main.f90 inputlist.f90 dfault.f90 inital.f90 setup.f90 output.f90 wrdump.f90 endrun.f90 numinc.f90 energy.f90 elapsed_time.f90 Modules.o 
#Modules.o:                     Modules.f90 
#ae_profiles.o:                 ae_profiles.f90 eqsplns.f90 functions.f90 Modules.o
#b2lx.o:                        b2lx.f90 Modules.o
#block.o:                       block.f90 Modules.o
#clgam.o:                       clgam.f90 Modules.o
#cnvt.o:                        cnvt.f90 Modules.o
#decbt.o:                       decbt.f90 Modules.o
#derivatives.o:                 derivatives.f90 Modules.o
#dfault.o:                      dfault.f90 Modules.o
#dlsq.o:                        dlsq.f90 eqdy_dyeq.f90 mult.f90 derivatives.f90 Modules.o
#dlstar.o:                      dlstar.f90 eqdy_dyeq.f90 mult.f90 derivatives.f90 Modules.o
#eigensolver_tools.o:           eigensolver.f90 Modules.o
#elapsed_time.o:                elapsed_time.f90 Modules.o
#endrun.o:                      endrun.f90 eqdy_dyeq.f90 mult.f90 derivatives.f90 decbt.f90 solbt.f90 Modules.o
#energy.o:                      energy.f90 eqdy_dyeq.f90 mult.f90 derivatives.f90 quadq.f90 Modules.o
#eqdy_dyeq.o:                   eqdy_dyeq.f90 Modules.o
#eqsplns.o:                     eqsplns.f90 eqsplns.f90 Modules.o
#etachi.o:                      etachi.f90 eqsplns.f90 functions.f90 Modules.o
#fitter.o:                      fitter.f90 eqsplns.f90 Modules.o
#functions.o:                   functions.f90 Modules.o
#grid.o:                        grid.f90 Modules.o
#inputlist.o:                   inputlist.f90 Modules.o
#Laundamp_tools.o:              Laundamp_tools.f90 Modules.o
#lincheck.o:                    lincheck.f90 clgam.f90 om.f90 b2lx.f90 derivatives.f90 Modules.o
#linstart.o:                    linstart.f90 cnvt.f90 clgam.f90 om.f90 block.f90 eigensolver.f90 derivatives.f90 decbt.f90 dlstar.f90 dlsq.f90 Modules.o
#linstep.o:                     linstep.f90 clgam.f90 b2lx.f90 solbt.f90 cnvt.f90 Modules.o
#mmlims.o:                      mmlims.f90 Modules.o
#mult.o:                        mult.f90 Modules.o
#numinc.o:                      numinc.f90 Modules.o
#om.o:                          om.f90 block.f90 Modules.o
#output.o:                      output.f90 eqdy_dyeq.f90 derivatives.f90 Modules.o
#pert.o:                        pert.f90 functions.f90 Modules.o
#quadq.o:                       quadq.f90 Modules.o
#rddump.o:                      rddump.f90 Modules.o
#seteq.o:                       seteq.f90 derivatives.f90 eqsplns.f90 quadq.f90 Modules.o
#setmod.o:                      setmod.f90 inputlist.f90 Modules.o
#setup.o:                       setup.f90 setmod.f90 mmlims.f90 grid.f90 seteq.f90 etachi.f90 pert.f90 Modules.o
#solbt.o:                       solbt.f90 Modules.o
#vmec.o:                        vmec.f90 mmlims.f90 mult.f90 derivatives.f90 eqsplns.f90 fitter.f90 Laundamp_tools.f90 Modules.o
#wrdump.o:                      wrdump.f90 Modules.o
#-----------------------------------------------------------------------                                                                                        
