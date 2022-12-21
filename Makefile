
CXX := $(shell root-config --cxx)
CXXFLAGS := $(shell root-config --cflags) -fPIC
LDFLAGS := $(shell root-config --glibs)

VERSION := $(shell root-config --version | cut -d'.' -f1)
ifeq ($(VERSION),6)
  DICTEXE = rootcling
  MINOR := $(shell root-config --version | cut -d'.' -f2 | cut -d'/' -f1)
  LT20 := $(shell [ $(MINOR) -lt 20 ] && echo true)
  ifeq ($(LT20),true)
    DICTFLAGS = -c -p
  endif
else
  DICTEXE = rootcint
  DICTFLAGS = -c -p
endif

PACKAGE = OscProb
HEADERS := $(filter-out $(CURDIR)/LinkDef.h, $(wildcard $(CURDIR)/*.h))
SOURCES := $(wildcard $(CURDIR)/*.cxx)
TARGET = lib$(PACKAGE)
TARGET_LIB = $(CURDIR)/$(TARGET).so
DICTIONARY = $(CURDIR)/tmp/$(TARGET).cxx

# GSL library
GSL_INCS := $(shell gsl-config --cflags)
GSL_LIBS := $(shell gsl-config --libs)

#Eigen library
Eigen_INCS = ${CURDIR}/eigen/Eigen

override CXXFLAGS += $(GSL_INCS) -I$(Eigen_INCS) -I$(CURDIR)/utils
override GSL_INCS += -I$(Eigen_INCS) -I${CURDIR}/utils
override LDFLAGS  += $(GSL_LIBS)

# the sets of directories to do various things in
MATRIX = $(CURDIR)/MatrixDecomp/libMatrixDecomp.so
SUBDIRS = MatrixDecomp
BUILDDIRS = $(SUBDIRS:%=build-%)
CLEANDIRS = $(SUBDIRS:%=clean-%)

PREMDIR = $(CURDIR)/PremTables
PREMFILE = $(PREMDIR)/prem_default.txt
PREMINC = $(CURDIR)/prem_default.hpp

all: $(Eigen_INCS) $(BUILDDIRS) $(PREMINC) $(TARGET_LIB)

$(TARGET_LIB): $(DICTIONARY) $(SOURCES) $(MATRIX)
	@echo "  Building $(PACKAGE)..."
	@$(CXX) $(CXXFLAGS) -O3 -shared -o$@ $^ $(LDFLAGS)

$(DICTIONARY): $(HEADERS) LinkDef.h
	@echo "  Making dictionary for $(PACKAGE)..."
	@mkdir -p tmp
	@$(DICTEXE) -f $@ $(DICTFLAGS) $(GSL_INCS) $^
	@if [ -e $(DICTIONARY:%.cxx=%_rdict.pcm) ] ; then mv -f $(DICTIONARY:%.cxx=%_rdict.pcm) $(CURDIR)/ ; fi # ROOT 6

$(Eigen_INCS):
	@echo "  Eigen not found. Trying to initialize submodule..."
	@echo "  git submodule update --init"
	@git submodule update --init

$(BUILDDIRS):
	@exec $(MAKE) -s -C $(@:build-%=%)

$(CLEANDIRS):
	@exec $(MAKE) -s -C $(@:clean-%=%) clean

$(PREMINC): $(PREMDIR) $(PREMFILE)
	@echo "#include <string>" > $@
	@echo "const std::string PREM_DIR = \"$(PREMDIR)\";" >> $@
	@echo "const std::string PREM_DEFAULT = \"$(PREMFILE)\";" >> $@

test: $(TARGET_LIB)
	@cd test && root -l -b -q ../tutorial/LoadOscProb.C ../tutorial/SetNiceStyle.C TestMethods.C
	@cd test && root -l -b -q ../tutorial/LoadOscProb.C ../tutorial/SetNiceStyle.C CheckProbs.C

clean: $(CLEANDIRS)
	@echo "  Cleaning $(PACKAGE)..."
	@rm -f $(TARGET_LIB)
	@rm -f $(TARGET_LIB:%.so=%_rdict.pcm) # ROOT 6
	@rm -f $(PREMINC)
	@rm -f $(DICTIONARY)
	@rm -f $(DICTIONARY:%.cxx=%.h) # ROOT 5
	@rm -rf tmp

.PHONY: $(SUBDIRS) $(BUILDDIRS) $(CLEANDIRS) test clean all
