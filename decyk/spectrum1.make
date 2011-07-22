# Makefile for spectrum1 and vspectrum1

GOBJS = nullgks1.o nullgks2.o gksnull.o
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
#GOBJS = libgks1.o libgks2.o libmcX.o
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
#GOBJS = libgks1.o libgks2.o libygl.o
#LIBS =  -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# Postcript printer
#GOBJS = libgks1.o libgks2.o libpsp.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs
# No graphics
#LIBS = -lSystemStubs

# Makefile gfortran compiler with MacOS X

#FC90 = gfortran
#FC77 = gfortran
#CC = gcc

#OPTS90 = -O3 -fdefault-real-8
#OPTS77 = -O3 -fdefault-real-8
#CCOPTS = -O
#MOPTS = -fno-automatic
#MBOPTS = -fno-automatic
#LOPTS =
#LEGACY =

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libgks1.o libgks2.o libygl.o
#LIBS = $(CARBON) -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# Postcript printer
#GOBJS = libgks1.o libgks2.o libpsp.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs
# TIFF graphics
#GOBJS = libgks1.o libgks2.o librstr.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs
# ncar graphics (single precision only) from http://ngwww.ucar.edu/
#GOBJS = libgks1.o libgks2.o ncarstub.o
#LIBS = $(CARBON) -lSystemStubs -L/$(NCARG_ROOT)/lib -lncarg_gks -lncarg_c \
#-L/usr/X11R6/lib -lX11
# Tektronix graphics
#GOBJS = libgks1.o libgks2.o libplt10.o plot10.o libt1.o libloc1.o
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
#GOBJS = libgks1.o libygl.o
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

FC90 = gfortran
FC77 = gfortran
CC = gcc

OPTS90 = -O3 -fdefault-real-8
OPTS77 = -O3 -fdefault-real-8
CCOPTS = -O
MOPTS = -fno-automatic
MBOPTS = -fno-automatic
LOPTS = -lpthread
LEGACY =

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
GOBJS = libgks1.o libgks2.o libygl.o
#LIBS = -L/usr/X11R6/lib -lYgl -lX11
LIBS = -L/usr/lib64 -lYgl -L/usr/X11R6/lib -lX11
# Postcript printer
#GOBJS = libgks1.o libgks2.o libpsp.o libloc1.o
#LIBS =
# TIFF graphics
#GOBJS = libgks1.o libgks2.o librstr.o libloc1.o
#LIBS =
# Tektronix graphics
#GOBJS = libgks1.o libgks2.o libplt10.o plot10.o libt1.o libloc1.o
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

OBJS = globals.o input1mod.o init1mod.o fft1mod.o field1mod.o diag1mod.o \
init1lib.o fft1lib.o field1lib.o diag1lib.o libanls1.o nlibpars.o

# Linkage rule

all : spectrum1 vspectrum1

spectrum1 : spectrum1.o $(OBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o spectrum1 spectrum1.o $(OBJS) \
        $(GOBJS) $(LIBS)

vspectrum1 : vspectrum1.o $(OBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o vspectrum1 vspectrum1.o $(OBJS) \
        $(GOBJS) $(LIBS)

# Compilation rules

gksnull.o : gksnull.f
	$(FC77) $(OPTS77) -c gksnull.f

libmcX.o : libmcX.f
	$(FC77) $(OPTS77) -c libmcX.f

libygl.o : libygl.f
	$(FC77) $(OPTS77) -c libygl.f

libpsp.o : libpsp.f
	$(FC77) $(OPTS77) -c libpsp.f

ncarstub.o : ncarstub.f
	$(FC77) $(OPTS77) -c ncarstub.f

libplt10.o : libplt10.f
	$(FC77) $(OPTS77) -c libplt10.f

plot10.o : plot10.f
	$(FC77) $(OPTS77) -c plot10.f

libt1.o : libt1.f
	$(FC77) $(OPTS77) -c libt1.f

libloc1.o : libloc1.f
	$(FC77) $(OPTS77) -c libloc1.f

libgks1.o : libgks1.f
	$(FC77) $(OPTS77) -c libgks1.f

libgks2.o : libgks2.f
	$(FC77) $(OPTS77) -c libgks2.f

nullgks1.o : nullgks1.f
	$(FC77) $(OPTS77) -c nullgks1.f

nullgks2.o : nullgks2.f
	$(FC77) $(OPTS77) -c nullgks2.f

init1lib.o : init1lib.f
	$(FC77) $(OPTS77) -c init1lib.f

fft1lib.o : fft1lib.f
	$(FC77) $(OPTS77) -c fft1lib.f

field1lib.o : field1lib.f
	$(FC77) $(OPTS77) -c field1lib.f

diag1lib.o : diag1lib.f
	$(FC77) $(OPTS77) -c diag1lib.f

libanls1.o : libanls1.f
	$(FC77) $(OPTS77) -c libanls1.f

nlibpars.o : nlibpars.f
	$(FC77) $(OPTS77) -c nlibpars.f

globals.o : globals.f
	$(FC90) $(OPTS90) -c globals.f

input1mod.o : input1mod.f globals.o
	$(FC90) $(OPTS90) -c input1mod.f

init1mod.o : init1mod.f input1mod.o
	$(FC90) $(OPTS90) -c init1mod.f

fft1mod.o : fft1mod.f diag1mod.o
	$(FC90) $(OPTS90) -c fft1mod.f

field1mod.o : field1mod.f globals.o
	$(FC90) $(OPTS90) $(LEGACY) -c field1mod.f

diag1mod.o : diag1mod.f init1mod.o
	$(FC90) $(OPTS90) -c diag1mod.f

spectrum1.o : spectrum1.f fft1mod.o field1mod.o diag1mod.o
	$(FC90) $(OPTS90) $(LEGACY) $(MBOPTS) -c spectrum1.f

vspectrum1.o : vspectrum1.f fft1mod.o field1mod.o diag1mod.o
	$(FC90) $(OPTS90) $(LEGACY) $(MBOPTS) -c vspectrum1.f

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f spectrum1 vspectrum1
