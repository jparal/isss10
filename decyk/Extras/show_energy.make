# Makefile for spectrum2 and vspectrum2

GOBJS = nullgks1.o
CARBON = /System/Library/Frameworks/Carbon.framework/Carbon

# Makefile Absoft compiler with MacOS X

#FC90 = f90
#FC77 = f77
#CC = gcc

#OPTS90 = -O3 -N113 -YEXT_NAMES=LCS -YEXT_SFX=_
#OPTS77 = -O -N113 -f -N15
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s -N11
#LOPTS = -plainappl
#LEGACY =

# Mac graphics
#GOBJS = libgks1.o libmcX.o
#LIBS = $(CARBON)
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libgks1.o libgks2.o libygl.o
#LIBS = -L/usr/X11R6/lib -lYgl -lX11
# Postcript printer
#GOBJS = libgks1.o libgks2.o libpsp.o libloc1.o
#LIBS =
# No graphics
#LIBS =

# Makefile Nag compiler with MacOS X

#FC90 = f95
#FC77 = f95
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY = -dusty

# No graphics
#LIBS =

# Makefile IBM compiler with MacOS X

#FC90 = xlf90
#FC77 = xlf
#CC = gcc

#OPTS90 = -O3 -qautodbl=dbl4 -qarch=g5
#OPTS77 = -O3 -qautodbl=dbl4 -qarch=g5 -qfixed
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#LEGACY =

# No graphics
#LIBS =

# Makefile g95 compiler with MacOS X

#FC90 = g95
#FC77 = g95
#CC = gcc

#OPTS90 = -O3 -r8 -fno-second-underscore
#OPTS77 = -O3 -r8 -fno-second-underscore
#CCOPTS = -O
#MOPTS = -fstatic
#MBOPTS = -fstatic
#LOPTS =
#LEGACY =

# X11 graphics
#GOBJS = libgks1.o libygl.o
#LIBS =  -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# Postcript printer
#GOBJS = libgks1.o libgks2.o libpsp.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs
# No graphics
#LIBS = -lSystemStubs

# Makefile gfortran compiler with MacOS X

FC90 = gfortran
FC77 = gfortran
CC = gcc

OPTS90 = -O3 -fdefault-real-8
OPTS77 = -O3 -fdefault-real-8
CCOPTS = -O
MOPTS = -fno-automatic
MBOPTS = -fno-automatic
LOPTS =
LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
GOBJS = libgks1.o libygl.o
LIBS = $(CARBON) -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# Postcript printer
#GOBJS = libgks1.o libpsp.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs
# TIFF graphics
#GOBJS = libgks1.o librstr.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs
# ncar graphics (single precision only) from http://ngwww.ucar.edu/
#GOBJS = libgks1.o ncarstub.o
#LIBS = $(CARBON) -lSystemStubs -L/$(NCARG_ROOT)/lib -lncarg_gks -lncarg_c \
#-L/usr/X11R6/lib -lX11
# Tektronix graphics
#GOBJS = libgks1.o libplt10.o plot10.o libt1.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs 
# No graphics
#LIBS = $(CARBON) -lSystemStubs

# Makefile Intel compiler with MacOS X

#FC90 = ifort
#FC77 = ifort
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS =
#LEGACY =

# X11 graphics
#GOBJS = dlibgks2.o libgks2.o libgks1.o libygl.o
#LIBS = -L/usr/X11R6/lib -lYgl -lX11
# Postcript printer
#GOBJS = libgks1.o libgks2.o libpsp.o libloc1.o
#LIBS = $(CARBON)
# No graphics
#LIBS =

# Makefile Intel compiler with Linux

#FC90 = ifort
#FC77 = ifort
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS =
#LEGACY =

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libgks1.o libgks2.o libygl.o
#LIBS = -L/usr/X11R6/lib -lYgl -lX11
#LIBS = -L/usr/lib64 -lYgl -L/usr/X11R6/lib -lX11
# Postcript printer
#GOBJS = libgks1.o libgks2.o libpsp.o libloc1.o
#LIBS =
# No graphics
#LIBS =

# Makefile gfortran compiler with Linux

#FC90 = gfortran
#FC77 = gfortran
#CC = gcc

#OPTS90 = -O3 -fdefault-real-8
#OPTS77 = -O3 -fdefault-real-8
#CCOPTS = -O
#MOPTS = -fno-automatic
#MBOPTS = -fno-automatic
#LOPTS = -lpthread
#LEGACY =

#MPOBJS = MacMPf77.o LnxMP.o LPProcessors.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libgks1.o libygl.o
#LIBS = -L/usr/X11R6/lib -lYgl -lX11
#LIBS = -L/usr/lib64 -lYgl -L/usr/X11R6/lib -lX11
# Postcript printer
#GOBJS = libgks1.o libpsp.o libloc1.o
#LIBS =
# TIFF graphics
#GOBJS = libgks1.o librstr.o libloc1.o
#LIBS =
# Tektronix graphics
#GOBJS = libgks1.o libplt10.o plot10.o libt1.o libloc1.o
#LIBS = 
# No graphics
#LIBS =

# Makefile Pathscale compiler with AMD Opteron and Linux

#FC90 = pathf90
#FC77 = pathf90
#CC = gcc

#OPTS90 = -O3 -r8 -static
#OPTS77 = -O3 -r8 -static
#CCOPTS = -O
#MOPTS = -static-data
#MBOPTS = -static-data
#LOPTS =
#LEGACY =

# No graphics
#LIBS =

#

OBJS = errors_class.o graf1d_class.o

# Linkage rule

all : show_es_energy.out show_em_energy.out show_em_momentum.out

show_es_energy.out : show_es_energy.o $(OBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o show_es_energy.out show_es_energy.o \
        $(OBJS) $(GOBJS) $(LIBS)
	
show_em_energy.out : show_em_energy.o $(OBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o show_em_energy.out show_em_energy.o \
        $(OBJS) $(GOBJS) $(LIBS)
	
show_em_momentum.out : show_em_momentum.o $(OBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o show_em_momentum.out show_em_momentum.o \
        $(OBJS) $(GOBJS) $(LIBS)

# Compilation rules

libmcX.o : libmcX.f
	$(FC77) $(OPTS77) -c libmcX.f

libygl.o : libygl.f
	$(FC77) $(OPTS77) -c libygl.f
	
libgks1.o : libgks1.f
	$(FC77) $(OPTS77) -c libgks1.f

nullgks1.o : nullgks1.f
	$(FC77) $(OPTS77) -c nullgks1.f

errors_class.o : errors_class.f
	$(FC90) $(OPTS90) -c errors_class.f

graf1d_class.o : graf1d_class.f errors_class.o
	$(FC90) $(OPTS90) -c graf1d_class.f


show_es_energy.o : show_es_energy.f graf1d_class.o
	$(FC90) $(OPTS90) $(LEGACY) $(MBOPTS) -c show_es_energy.f

show_em_energy.o : show_em_energy.f graf1d_class.o
	$(FC90) $(OPTS90) $(LEGACY) $(MBOPTS) -c show_em_energy.f

show_em_momentum.o : show_em_momentum.f graf1d_class.o
	$(FC90) $(OPTS90) $(LEGACY) $(MBOPTS) -c show_em_momentum.f

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f *.out
