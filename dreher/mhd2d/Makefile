SHELL=/bin/sh

SRC := $(wildcard *.cc)
OBJ := $(SRC:%.cc=%.o)
DEP := $(SRC:%.cc=%.d)
EXE := mhd2D

CXX := g++
CXXFLAGS := -O2 -Wall

auto:  $(EXE)

info:
	echo $(SRC)

new:  clean auto

clean: 
	rm -f $(OBJ) $(DEP) $(EXE)

$(OBJ): %.o: %.cc
	 $(CXX) $(CXXFLAGS) -c $<

$(DEP): %.d: %.cc
	@$(SHELL) -ec 'g++ -MM  $< | sed '\''s/\($*\.o\)/\1 $@/g'\'' > $@'

$(EXE): $(OBJ)
	$(CXX) -o $(EXE) $(OBJ)  $(LINKOPTS) $(LIB_MPI)

#########################################################
# include dependency files
#########################################################
ifdef DEP
-include $(DEP)
endif
