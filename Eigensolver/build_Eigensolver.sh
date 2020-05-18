#!/bin/bash
#-----------------------------------------------------------------------
# Makefile for FAR3d Eigensolver executable
#-----------------------------------------------------------------------
#This has been modified from Jacobos distribution for the Mac with ifort compiler
echo Welcome to FAR3d Eigensolver

#These are generic flags and commands.
Comp="ifort" 
Flag=" -fixed"
MKLINCLUDE="-I${MKLROOT}/include"
MKLLINK="${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_core.a ${MKLROOT}/lib/libmkl_sequential.a -lpthread -lm"

#Comp="gfortran" 
#Flag=" -ffixed-form"
#Flag=" -ffixed-form -Wall -fcheck=all -g -fbacktrace -Wno-tabs"
#Flag=" -ffixed-form -Wall -g -frange-check -fbounds-check -Wcharacter-truncation -fcheck=all -fbacktrace -ffpe-trap=invalid -ftrapv -Waliasing -Wampersand -Wconversion -Wsurprising -Wc-binding-type -Wintrinsics-std -Wintrinsic-shadow -Wline-truncation -Wtarget-lifetime -Wno-tabs"

Exc=" xEigen" 
Opt=" -O2"

echo The Eigensolver executable is been compiled

#Libraries:
OBJS=" LIB_JDQZ/error.f LIB_JDQZ/jdqz.f LIB_JDQZ/jdqzmv.f LIB_JDQZ/makemm.f LIB_JDQZ/mkqkz.f LIB_JDQZ/myexc.f" 
OBJS2=" LIB_JDQZ/psolve.f LIB_JDQZ/qzsort.f LIB_JDQZ/select.f LIB_JDQZ/zcgstab.f LIB_JDQZ/zgmres.f"
OBJS3=" LIB_JDQZ/zmgs.f LIB_JDQZ/zones.f LIB_JDQZ/zxpay.f LIB_JDQZ/zzeros.f"

#Main eigensolver files:
OBJS4=" grandom_mod.f90 symtrd_mod.f90 myjdqz_mod_TAEFL_cmplx.f90 main_jdqz_TAEFL_cmplx.f90"

#lapack library:
#OBJS5=" -I. -lblas -llapack"
OBJS5=" -I. $MKLINCLUDE $MKLLINK"

$Comp -o $Exc $Opt $Flag $OBJS $OBJS2 $OBJS3 $OBJS4 $OBJS5

echo The FAR3D Eigensolver executable has been created
