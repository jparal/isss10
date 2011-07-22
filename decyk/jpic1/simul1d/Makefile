#Makefile for 1D PIC codes in new_beps1.source
# Java Package simul1d for functions in simul1mod.f

JCOPTS = -g -deprecation -source 5

all: globals.class push1d.class fft1d.class field1d.class init1d.class \
     input1d.class diag1d.class bpush1d.class dpush1d.class

globals.class : globals.java
	javac $(JCOPTS) globals.java

push1d.class : push1d.java
	javac $(JCOPTS) push1d.java

fft1d.class : fft1d.java
	javac $(JCOPTS) fft1d.java

field1d.class : field1d.java
	javac $(JCOPTS) field1d.java

init1d.class : init1d.java
	javac $(JCOPTS) init1d.java

input1d.class : input1d.java globals.class
	javac $(JCOPTS) input1d.java

diag1d.class : diag1d.java
	javac $(JCOPTS) diag1d.java

bpush1d.class : bpush1d.java
	javac $(JCOPTS) bpush1d.java

dpush1d.class : dpush1d.java
	javac $(JCOPTS) dpush1d.java

clean :
	rm -f *.class

