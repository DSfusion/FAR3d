FC     = ifort
FFLAGS = -fixed
FOPT   = -O2

OBJ    = error.o jdqz.o jdqzmv.o makemm.o mkqkz.o myexc.o \
         psolve.o qzsort.o select.o zcgstab.o zgmres.o \
         zmgs.o zones.o zxpay.o zzeros.o

%.o: %.f
	$(FC) $(FFLAGS) $(FOPT) -c $<

libjdqz.a: $(OBJ)
	ar -r $@ $(OBJ)

clean:
	rm -f *.o libjdqz.a
