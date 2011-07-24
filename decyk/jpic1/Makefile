#Makefile for 1D PIC codes in new_beps1.source
# Java Interface

# Mac OS X
#JHEADERS1 = /System/Library/Frameworks/JavaVM.framework/Headers
#JHEADERS2 =
# OpenJDK Linux
JHEADERS1 = /usr/lib/jvm/java-1.6.0-openjdk-1.6.0.0.x86_64/include
JHEADERS2 = -I /usr/lib/jvm/java-1.6.0-openjdk-1.6.0.0.x86_64/include/linux

CC = gcc
CCOPTS = -O -fPIC
JCOPTS = -g -deprecation -source 5 -classpath .

COBJS = push1lib_f.o fft1lib_f.o field1lib_f.o init1lib_f.o diag1lib_f.o \
bpush1lib_f.o dpush1lib_f.o dtimer.o

JOBJS = push1d_j.o fft1d_j.o field1d_j.o init1d_j.o diag1d_j.o \
bpush1d_j.o dpush1d_j.o

HOBJS = simul1d_push1d.h simul1d_fft1d.h simul1d_field1d.h simul1d_init1d.h \
simul1d_diag1d.h simul1d_bpush1d.h simul1d_dpush1d.h

all: libsimul1d.a

libsimul1d.a : $(JOBJS) $(COBJS)
	ar rc libsimul1d.a $(JOBJS) $(COBJS)

push1lib_f.o : push1lib_f.c
	$(CC) $(CCOPTS) -c push1lib_f.c

fft1lib_f.o : fft1lib_f.c
	$(CC) $(CCOPTS) -c fft1lib_f.c

field1lib_f.o : field1lib_f.c
	$(CC) $(CCOPTS) -c field1lib_f.c

init1lib_f.o : init1lib_f.c
	$(CC) $(CCOPTS) -c init1lib_f.c

diag1lib_f.o : diag1lib_f.c
	$(CC) $(CCOPTS) -c diag1lib_f.c
	
bpush1lib_f.o : bpush1lib_f.c
	$(CC) $(CCOPTS) -c bpush1lib_f.c

dpush1lib_f.o : dpush1lib_f.c
	$(CC) $(CCOPTS) -c dpush1lib_f.c

dtimer.o : dtimer.c
	$(CC) $(CCOPTS) -c dtimer.c

push1d_j.o : push1d_j.c simul1d_push1d.h
	$(CC) $(CCOPTS) -I $(JHEADERS1) $(JHEADERS2) -c push1d_j.c

fft1d_j.o : fft1d_j.c simul1d_fft1d.h
	$(CC) $(CCOPTS) -I $(JHEADERS1) $(JHEADERS2) -c fft1d_j.c

field1d_j.o : field1d_j.c simul1d_field1d.h
	$(CC) $(CCOPTS) -I $(JHEADERS1) $(JHEADERS2) -c field1d_j.c

init1d_j.o : init1d_j.c simul1d_init1d.h
	$(CC) $(CCOPTS) -I $(JHEADERS1) $(JHEADERS2) -c init1d_j.c

diag1d_j.o : diag1d_j.c simul1d_diag1d.h
	$(CC) $(CCOPTS) -I $(JHEADERS1) $(JHEADERS2) -c diag1d_j.c

bpush1d_j.o : bpush1d_j.c simul1d_bpush1d.h
	$(CC) $(CCOPTS) -I $(JHEADERS1) $(JHEADERS2) -c bpush1d_j.c

dpush1d_j.o : dpush1d_j.c simul1d_dpush1d.h
	$(CC) $(CCOPTS) -I $(JHEADERS1) $(JHEADERS2) -c dpush1d_j.c

simul1d_push1d.h : jlibs.class
	javah -jni simul1d.push1d

simul1d_fft1d.h : jlibs.class
	javah -jni simul1d.fft1d

simul1d_field1d.h : jlibs.class
	javah -jni simul1d.field1d

simul1d_init1d.h : jlibs.class
	javah -jni simul1d.init1d

simul1d_diag1d.h : jlibs.class
	javah -jni simul1d.diag1d

simul1d_bpush1d.h : jlibs.class
	javah -jni simul1d.bpush1d

simul1d_dpush1d.h : jlibs.class
	javah -jni simul1d.dpush1d

jlibs.class : jlibs.java simul1d
	javac $(JCOPTS) jlibs.java

beps1.class : beps1.java
	javac $(JCOPTS) beps1.java

bbeps1.class : bbeps1.java
	javac $(JCOPTS) bbeps1.java

dbeps1.class : dbeps1.java
	javac $(JCOPTS) dbeps1.java

jclean :
	make -C simul1d clean
	
clean : jclean
	rm -f *.class *.o *.a $(HOBJS)

