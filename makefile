##
## Makefile for KtJet examples.
##
## B. Waugh 2005-07-22
##

# Default executble to compile
EXECUTABLE = main
# Location of KtJet installation
#KTDIR = /data/zenith234b/kuprash/programs/ktjet_files
#KTDIR = /afs/desy.de/group/zeus/pool/kuprash/programs/KtJet-1.08
KTDIR = /afs/desy.de/user/g/glusheno/programs/KtJet-1.08
#KTDIR = /usr/local
# Location of CLHEP installation
#CLHEPDIR = /data/zenith234b/kuprash/programs/clhep_files
#CLHEPDIR = /afs/cern.ch/sw/lcg/external/clhep/2.0.4.0/slc4_amd64_gcc34
CLHEPDIR = /cvmfs/sft.cern.ch/lcg/external/clhep/2.0.4.5/x86_64-slc5-gcc43-opt
#CLHEPDIR = /usr/local
# Following line may be needed if using version 1 of CLHEP
#CLHEP1FLAG = -DCLHEP1
# Uncomment if library compiled with double precision
#KTFLAGS = -DKTDOUBLEPRECISION

# Programs and flags
CXX          = g++ 
CXXFLAGS     = -O0 -ansi -std=gnu++0x
RM           = rm -f 

# Assign or Add variables
CXXFLAGS    +=  -I$(KTDIR)/include -I$(CLHEPDIR)/include -Iinc $(shell root-config --cflags)
CXXFLAGS    += $(CLHEP1FLAG) $(KTFLAGS)
LIBS         = -L$(KTDIR)/lib -lKtJet -L$(CLHEPDIR)/lib -lCLHEP $(shell root-config --libs)
GARBAGE      = *.o src/*.o
SOURCE_FILES = $(shell ls -1 src/*.c)
SOURCE_FILESF:= $(shell ls -1 src/*.f)

OBJ          = $(patsubst %.c,%.o,$(SOURCE_FILES))
OBJ         += $(patsubst %.f,%.o,$(SOURCE_FILESF))
INC          = -I./inc


#vpath       %.h inc
.PHONY: clean
.PRECIOUS: %.o

default: $(EXECUTABLE) 

% : %.o $(OBJ) makefile
	$(CXX) -o $@ $@.o $(OBJ) $(CXXFLAGS) $(LIBS)  -lgfortran -B /lib/ssa

%.o : %.cxx 
	$(CXX) $(CXXFLAGS)  -c $<

src/selector.o : src/selector.c inc/selector.h inc/constants.h 
	$(CXX) $(CXXFLAGS) -o $@ -c src/selector.c

src/selectorBegin.o : src/selectorBegin.c inc/selector.h inc/constants.h
	$(CXX) $(CXXFLAGS) -o $@ -c src/selectorBegin.c

src/selectorTerminate.o : src/selectorTerminate.c inc/selector.h inc/constants.h
	$(CXX) $(CXXFLAGS) -o $@ -c src/selectorTerminate.c

%.o : %.c
	$(CXX) $(CXXFLAGS) $(INC) -o $@ -c $<

.f.o:
	g++ -std=gnu++0x -o $@ -c $<


#src/runinfo.o : src/runinfo.c inc/runinfo.h
#	$(CXX) $(CXXFLAGS) -o $@ -c  src/runinfo.c

#src/fill_chain.o : src/fill_chain.c inc/*Names.h
#	$(CXX) $(CXXFLAGS) -o $@ -c  src/fill_chain.c
#
#src/hist.o : src/hist.c inc/hist.h
#	$(CXX) $(CXXFLAGS) -o $@ -c  src/hist.c

clean :
	$(RM) $(GARBAGE) main 
