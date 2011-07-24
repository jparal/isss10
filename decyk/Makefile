#Makefile for 1D PIC codes in new_beps1.source

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
#NOSAVE =
#LEGACY =

#MPOBJS = LnxMP.o

# Mac graphics
#GOBJS = libgks1.o libmcX.o
#LIBS = $(CARBON)
# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libgks1.o libygl.o
#LIBS = $(CARBON)  -L/usr/X11R6/lib -lYgl -lX11
# No graphics
#LIBS = $(CARBON)

# Makefile Nag compiler with MacOS X

#FC90 = f95
#FC77 = f95
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS = -framework carbon 
#NOSAVE =
#LEGACY = -dusty

#MPOBJS = LnxMP.o

# No graphics
#LIBS =

# Makefile IBM compiler with MacOS X

#FC90 = xlf90
#FC77 = xlf
#CC = gcc

#OPTS90 = -O3 -qautodbl=dbl4 -qarch=g5 -qextname
#OPTS77 = -O3 -qautodbl=dbl4 -qarch=g5 -qextname -qfixed
#CCOPTS = -O
#MOPTS = -s
#MBOPTS = -s
#LOPTS =
#NOSAVE = -qnosave
#LEGACY =

#MPOBJS = LnxMP.o

# No graphics
#LIBS = $(CARBON)

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
#NOSAVE =
#LEGACY =

#MPOBJS = LnxMP.o

# X11 graphics
#GOBJS = libgks1.o libygl.o
#LIBS = $(CARBON) -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
# Postcript printer
#GOBJS = libgks1.o libpsp.o libloc1.o
#LIBS = $(CARBON) -lSystemStubs
# No graphics
#LIBS = $(CARBON) -lSystemStubs

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
#NOSAVE = -frecursive
#LEGACY =

#MPOBJS = LnxMP.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libgks1.o libygl.o
#LIBS = $(CARBON) -lSystemStubs -L/usr/X11R6/lib -lYgl -lX11
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

# Makefile Intel compiler with Mac OS X

#FC90 = ifort
#FC77 = ifort
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS = -lpthread
#NOSAVE = -automatic
#LEGACY =

#MPOBJS = LnxMP.o

# X11 graphics
#GOBJS = libgks1.o libygl.o
#LIBS = $(CARBON) -L/usr/X11R6/lib -lYgl -lX11
# Postcript printer
#GOBJS = libgks1.o libpsp.o libloc1.o
#LIBS = $(CARBON)
# No graphics
#LIBS = $(CARBON)

#CUDA
#LOPTS = -lpthread -cxxlib
#LIBS = $(CARBON) -L/usr/local/cuda/lib -lcuda -lcudart

# Makefile Intel compiler with Linux

#FC90 = ifort
#FC77 = ifort
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -save
#MBOPTS = -save
#LOPTS = -lpthread
#NOSAVE = -automatic
#LEGACY =

#MPOBJS = LnxMP.o LPProcessors.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
#GOBJS = libgks1.o libygl.o
#LIBS = -L/usr/X11R6/lib -lYgl -lX11
#LIBS = -L/usr/lib64 -lYgl -L/usr/X11R6/lib -lX11
# Postcript printer
#GOBJS = libgks1.o libpsp.o libloc1.o
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
NOSAVE = -frecursive
LEGACY =

MPOBJS = LnxMP.o LPProcessors.o

# X11 graphics needs Ygl from http://www.thp.uni-duisburg.de/Ygl/
GOBJS = libgks1.o libygl.o
#LIBS = -L/usr/X11R6/lib -lYgl -lX11
LIBS = -L/usr/lib64 -lYgl -L/usr/X11R6/lib -lX11
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

#OPTS90 = -O3 -r8 -fno-second-underscore
#OPTS77 = -O3 -r8 -fno-second-underscore
#CCOPTS = -O
#MOPTS = -static-data
#MBOPTS = -static-data
#LOPTS = -lpthread
#NOSAVE =
#LEGACY =

#MPOBJS = LnxMP.o LPProcessors.o

# No graphics
#LIBS =

# Makefile PGI compiler with AMD Opteron and Linux

#FC90 = pgf90
#FC77 = pgf90
#CC = gcc

#OPTS90 = -O3 -r8
#OPTS77 = -O3 -r8
#CCOPTS = -O
#MOPTS = -Msave
#MBOPTS = -Msave
#LOPTS = -lpthread
#NOSAVE =
#LEGACY =

#MPOBJS = LnxMP.o LPProcessors.o

# No graphics
#LIBS =

#

GESOBJS = espush1mod.o mp0mod.o
MGESOBJS = mespush1mod.o mp0mod.o

ESOBJS = globals.o input1mod.o init1mod.o diag1mod.o bpush1mod.o push1mod.o \
fft1mod.o field1mod.o init1lib.o bpush1lib.o push1lib.o fft1lib.o field1lib.o \
diag1lib.o

MESOBJS = mpush1mod.o mpush1lib.o

EMOBJS = dpush1mod.o dpush1lib.o

DESOBJS = dfield1mod.o dfield1lib.o

NPOBJS = nullMP.o

# Linkage rule

all : new_beps1 new_bbeps1 new_dbeps1

threaded : new_mbeps1gl

bounded: new_d0_beps1 new_d0_bbeps1 new_d0_dbeps1

new_beps1 : new_beps1.o $(GESOBJS) $(ESOBJS) $(NPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_beps1 \
        new_beps1.o $(GESOBJS) $(ESOBJS) $(NPOBJS) $(GOBJS) $(LIBS)

new_bbeps1 : new_bbeps1.o $(ESOBJS) $(EMOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_bbeps1 \
        new_bbeps1.o $(ESOBJS) $(EMOBJS) $(GOBJS) $(LIBS)

new_dbeps1 : new_dbeps1.o $(ESOBJS) $(EMOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_dbeps1 \
        new_dbeps1.o $(ESOBJS) $(EMOBJS) $(GOBJS) $(LIBS)

new_beps1gl : new_beps1gl.o $(GESOBJS) $(ESOBJS) $(NPOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_beps1gl \
        new_beps1gl.o $(GESOBJS) $(ESOBJS) $(NPOBJS) $(GOBJS) $(LIBS)

new_d0_beps1 : new_d0_beps1.o $(GESOBJS) $(ESOBJS) $(DESOBJS) $(NPOBJS) \
                   $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_d0_beps1 \
        new_d0_beps1.o $(GESOBJS) $(ESOBJS) $(DESOBJS) $(NPOBJS) $(GOBJS) \
    $(LIBS)

new_d0_bbeps1 : new_d0_bbeps1.o $(ESOBJS) $(EMOBJS) $(DESOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_d0_bbeps1 \
        new_d0_bbeps1.o $(ESOBJS) $(EMOBJS) $(DESOBJS) $(GOBJS) $(LIBS)

new_d0_dbeps1 : new_d0_dbeps1.o $(ESOBJS) $(EMOBJS) $(DESOBJS) $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_d0_dbeps1 \
        new_d0_dbeps1.o $(ESOBJS) $(EMOBJS) $(DESOBJS) $(GOBJS) $(LIBS)

new_mbeps1gl : new_mbeps1gl.o $(MGESOBJS) $(ESOBJS) $(MESOBJS) $(MPOBJS) \
                   $(GOBJS)
	$(FC90) $(OPTS90) $(LOPTS) -o new_mbeps1gl \
        new_mbeps1gl.o $(MGESOBJS) $(ESOBJS) $(MESOBJS) $(MPOBJS) $(GOBJS) \
    $(LIBS)

# Compilation rules

LnxMP.o : LnxMP.c
	$(CC) $(CCOPTS) -c LnxMP.c

LPProcessors.o : LPProcessors.c
	$(CC) $(CCOPTS) -c LPProcessors.c

libmcX.o : libmcX.f
	$(FC77) $(OPTS77) -c libmcX.f

libygl.o : libygl.f
	$(FC77) $(OPTS77) -c libygl.f

libpsp.o : libpsp.f
	$(FC77) $(OPTS77) -c libpsp.f

librstr.o : librstr.f
	$(FC77) $(OPTS77) -c librstr.f

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

nullgks1.o : nullgks1.f
	$(FC77) $(OPTS77) -c nullgks1.f
	
nullMP.o : nullMP.f
	$(FC77) $(OPTS77) -c nullMP.f

init1lib.o : init1lib.f
	$(FC77) $(OPTS77) -c init1lib.f

push1lib.o : push1lib.f
	$(FC77) $(OPTS77) -c push1lib.f

mpush1lib.o : mpush1lib.f
	$(FC77) $(OPTS77) $(LEGACY) $(NOSAVE) -c mpush1lib.f

bpush1lib.o : bpush1lib.f
	$(FC77) $(OPTS77) -c bpush1lib.f

dpush1lib.o : dpush1lib.f
	$(FC77) $(OPTS77) -c dpush1lib.f

fft1lib.o : fft1lib.f
	$(FC77) $(OPTS77) -c fft1lib.f

field1lib.o : field1lib.f
	$(FC77) $(OPTS77) -c field1lib.f

dfield1lib.o : dfield1lib.f
	$(FC77) $(OPTS77) -c dfield1lib.f

diag1lib.o : diag1lib.f
	$(FC77) $(OPTS77) -c diag1lib.f

globals.o : globals.f
	$(FC90) $(OPTS90) -c globals.f

input1mod.o : input1mod.f globals.o
	$(FC90) $(OPTS90) -c input1mod.f

init1mod.o : init1mod.f input1mod.o
	$(FC90) $(OPTS90) -c init1mod.f

mpush1mod.o : mpush1mod.f push1mod.o mp0mod.o
	$(FC90) $(OPTS90) -c mpush1mod.f

push1mod.o : push1mod.f diag1mod.o
	$(FC90) $(OPTS90) -c push1mod.f

bpush1mod.o : bpush1mod.f diag1mod.o
	$(FC90) $(OPTS90) -c bpush1mod.f

dpush1mod.o : dpush1mod.f diag1mod.o
	$(FC90) $(OPTS90) -c dpush1mod.f

fft1mod.o : fft1mod.f diag1mod.o
	$(FC90) $(OPTS90) -c fft1mod.f

mespush1mod.o : mespush1mod.f mpush1mod.o
	$(FC90) $(OPTS90) -c mespush1mod.f

espush1mod.o : espush1mod.f push1mod.o
	$(FC90) $(OPTS90) -c espush1mod.f

field1mod.o : field1mod.f globals.o
	$(FC90) $(OPTS90) $(LEGACY) -c field1mod.f

dfield1mod.o : dfield1mod.f globals.o
	$(FC90) $(OPTS90) -c dfield1mod.f

diag1mod.o : diag1mod.f init1mod.o
	$(FC90) $(OPTS90) -c diag1mod.f

mp0mod.o : mp0mod.f
	$(FC90) $(OPTS90) -c mp0mod.f

new_beps1.o : new_beps1.f espush1mod.o fft1mod.o field1mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_beps1.f

new_dbeps1.o : new_dbeps1.f bpush1mod.o dpush1mod.o push1mod.o \
               fft1mod.o field1mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_dbeps1.f

new_bbeps1.o : new_bbeps1.f bpush1mod.o dpush1mod.o push1mod.o \
               fft1mod.o field1mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_bbeps1.f

new_beps1gl.o : new_beps1gl.f espush1mod.o fft1mod.o field1mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_beps1gl.f

new_d0_beps1.o : new_d0_beps1.f espush1mod.o fft1mod.o field1mod.o dfield1mod.o \
                 mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_d0_beps1.f

new_d0_dbeps1.o : new_d0_dbeps1.f bpush1mod.o dpush1mod.o push1mod.o \
               fft1mod.o field1mod.o dfield1mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_d0_dbeps1.f

new_d0_bbeps1.o : new_d0_bbeps1.f bpush1mod.o dpush1mod.o push1mod.o \
                  fft1mod.o field1mod.o dfield1mod.o
	$(FC90) $(OPTS90) $(MOPTS) -c new_d0_bbeps1.f

new_mbeps1gl.o : new_beps1gl.f mespush1mod.o fft1mod.o field1mod.o mp0mod.o
	$(FC90) $(OPTS90) $(MOPTS) -o new_mbeps1gl.o -c new_beps1gl.f

clean :
	rm -f *.o *.mod

clobber: clean
	rm -f new_beps1 new_bbeps1 new_dbeps1 new_mbeps1gl new_d0_beps1 \
	      new_d0_bbeps1 new_d0_dbeps1
