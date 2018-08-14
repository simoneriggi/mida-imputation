#
# General Makefile for the OfflineUser package
#
#

# Replace the wildcard expression with .cc file list if you do
# not want to compile all .cc files in this directory
#
TOPDIR= $(shell pwd)

BINDIR  = $(TOPDIR)/bin
LIBDIR  = $(TOPDIR)/lib
SRCDIR  = $(TOPDIR)/src
INCDIR  = $(TOPDIR)/include
OBJDIR  = $(TOPDIR)/obj

USER_SRCS = $(wildcard $(SRCDIR)/*.cc)
R_HOME = $(shell R RHOME)
HEADERS_DICT = $(INCDIR)/DataReader.h $(INCDIR)/Util.h $(INCDIR)/MathUtils.h  $(INCDIR)/MultipleImputation.h $(INCDIR)/MeanImputation.h $(INCDIR)/ListwiseDeletion.h $(INCDIR)/KMeansClustering.h $(INCDIR)/MNMixtureClustering.h


OBJS = $(USER_SRCS:.cc=.o)

## Get platform 32/64 bit
LBITS   = $(shell getconf LONG_BIT)

# Set executable a name
ROOT_CLASS_DICT= ClassDictionary.cc
SHARED_LIB = libMDImputation.so
MAIN_SRC= MDImputation.cc
EXE = MDImputation
#
#############################################################

## You should not need to change anything below this line ###

.PHONY: all depend clean


######################################
###  CPPFLAGS & CXXFLAGS  ############
######################################

CPPFLAGS = -I$(INCDIR)

ifeq ($(LBITS),64)
  # do 64 bit stuff here
	CPPFLAGS += -I/usr/include -pthread -m64
	CXXFLAGS = -std=c++11 -O2 -g -Wall -Wunused -Wuninitialized -Woverloaded-virtual -fPIC -pthread -m64 
	SYSLIBDIR = /usr/lib64
else
  # do 32 bit stuff here
	CPPFLAGS += -I/usr/include -pthread -m32
	CXXFLAGS = -std=c++11 -O2 -g -Wall -Wunused -Wuninitialized -Woverloaded-virtual -fPIC -pthread -m32 
	SYSLIBDIR = /usr/lib
endif

SOFLAGS = -fPIC -ggdb3 -Wall -shared

CPPFLAGS_ROOT= -I$(ROOTSYS)/include

CPPFLAGS_R = $(shell $(R_HOME)/bin/R CMD config --cppflags)
CPPFLAGS_R += $(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
CPPFLAGS_RINSIDE = $(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)

CPPFLAGS += $(CPPFLAGS_ROOT) 
CPPFLAGS += $(CPPFLAGS_R) $(CPPFLAGS_RINSIDE) 

###########################
###  LDFLAGS   ############
###########################
LDFLAGS_ROOT = $(shell root-config --libs) -L$(ROOTSYS)/lib -lSpectrum -lMathMore -lMinuit
LDFLAGS_SYSTEM = -L$(SYSLIBDIR) -lrt

LDFLAGS_R  = $(shell $(R_HOME)/bin/R CMD config --ldflags)
LDFLAGS_R += $(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
LDFLAGS_R += $(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)
LDFLAGS_R += $(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)
LDFLAGS_RINSIDE = $(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

LDFLAGS = $(LDFLAGS_SYSTEM)
LDFLAGS += $(LDFLAGS_ROOT) 
LDFLAGS += $(LDFLAGS_R) $(LDFLAGS_RINSIDE)
################################################################


all: GETOBJS $(SHARED_LIB) $(EXE) PUTOBJS


PRINTINFO: 
	@echo 'Compiling $(EXE) on a $(LBITS) bit machine' \

GETOBJS:
	@echo "Put object files again to $(SRCDIR) dir"
	- mv -f $(OBJDIR)/*.o $(SRCDIR) 2>/dev/null
	- mv -f $(LIBDIR)/$(SHARED_LIB) $(TOPDIR) 2>/dev/null
	- mv -f $(LIBDIR)/ClassDictionary_rdict.pcm $(TOPDIR) 2>/dev/null
	- mv -f $(BINDIR)/$(EXE) $(TOPDIR) 2>/dev/null

PUTOBJS:
	@echo "Moving object files to $(OBJDIR) dir"
	- mv -f $(SRCDIR)/*.o $(OBJDIR) 2>/dev/null
	- mv -f *.o $(OBJDIR) 2>/dev/null
	- mv -f $(SHARED_LIB) $(LIBDIR) 2>/dev/null
	- mv -f ClassDictionary_rdict.pcm $(LIBDIR) 2>/dev/null
	- mv -f $(EXE) $(BINDIR) 2>/dev/null

$(ROOT_CLASS_DICT): $(HEADERS_DICT) LinkDef.h
	@echo "Generating $@ ROOT dictionary class ..."
	rootcling -f $@ -c $(CXXFLAGS) $(CPPFLAGS) -p $^

$(SHARED_LIB): $(ROOT_CLASS_DICT) $(USER_SRCS)
	@echo "Generating $@ shared library ..."
	@$(CXX) $(CXXFLAGS) $(SOFLAGS) $(CPPFLAGS) $^ -o $@ $(LDFLAGS) 

$(EXE): $(MAIN_SRC) $(OBJS) ClassDictionary.o
	@echo "Building $(EXE) on a $(LBITS) bit machine"
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDFLAGS) -L$(TOPDIR) -lMDImputation



#############################################################
# gcc can generate the dependency list

depend: Make-depend

Make-depend: $(USER_SRCS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -MM $^ > $@


clean:
	- rm -f *.o $(OBJDIR)/*.o $(SRCDIR)/*.o
	- rm -f $(EXE) $(BINDIR)/$(EXE)
	- rm -f $(SHARED_LIB) $(LIBDIR)/$(SHARED_LIB)
	- rm -f ClassDictionary.cc ClassDictionary.h ClassDictionary_rdict.pcm $(LIBDIR)/ClassDictionary_rdict.pcm
	- rm -f core Make-depend

include Make-depend

