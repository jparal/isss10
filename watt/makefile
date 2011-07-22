F77 = gfortran
CC = cc
FFLAGS = -C -lc 
CFLAGS = -c 
DEBUG =

FOBJS1 = v1d1code.f90
FOBJS2 = v1d1simcode.f90

all:	v1d1.x v1d1sim.x 

v1d1.x:	$(FOBJS1) 
		$(F77) $(FOBJS1) -g -o v1d1.x

v1d1sim.x:      $(FOBJS2) 
		$(F77) $(FOBJS2) -g -o v1d1sim.x



