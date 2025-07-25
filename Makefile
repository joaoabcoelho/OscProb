
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
HEADERS := $(filter-out $(CURDIR)/inc/LinkDef.h, $(wildcard $(CURDIR)/inc/*.h))
SOURCES := $(wildcard $(CURDIR)/src/*.cxx)
TARGET = lib$(PACKAGE)
TARGET_LIB = $(CURDIR)/lib/$(TARGET).so
TARGET_PCM = $(CURDIR)/lib/$(TARGET)_rdict.pcm
DICTIONARY = $(CURDIR)/tmp/$(TARGET).cxx

#Eigen library
Eigen_INCS = ${CURDIR}/eigen

INCDIRS = -I$(CURDIR) -I$(CURDIR)/inc -I$(Eigen_INCS)

override CXXFLAGS += $(INCDIRS)

# the sets of directories to do various things in
MATRIX = $(CURDIR)/MatrixDecomp/libMatrixDecomp.so
SUBDIRS = MatrixDecomp
BUILDDIRS = $(SUBDIRS:%=build-%)
CLEANDIRS = $(SUBDIRS:%=clean-%)

PREMDIR = $(CURDIR)/PremTables
MODEL3DDIR = $(CURDIR)/EarthTables
PREMFILE = $(PREMDIR)/prem_default.txt
PREM3DFILE = $(MODEL3DDIR)/earth_binned_default.txt
PREMINC = $(CURDIR)/inc/prem_default.hpp

# Define list of submodules
SUBMODULES = $(shell grep path .gitmodules | sed 's/.*= //')
SUBMODULES := $(patsubst %,%/.git,$(SUBMODULES))

all: $(SUBMODULES) $(BUILDDIRS) $(PREMINC) $(TARGET_LIB)

$(TARGET_LIB): $(DICTIONARY) $(SOURCES) $(MATRIX)
	@echo "  Building $(PACKAGE)..."
	@mkdir -p lib
	@$(CXX) $(CXXFLAGS) -O3 -shared -o$@ $^ $(LDFLAGS)

$(DICTIONARY): $(HEADERS) inc/LinkDef.h
	@echo "  Making dictionary for $(PACKAGE)..."
	@mkdir -p tmp
	@$(DICTEXE) -f $@ $(DICTFLAGS) $(INCDIRS) $^
	@if [ -e $(DICTIONARY:%.cxx=%_rdict.pcm) ] ; then mkdir -p lib && mv $(DICTIONARY:%.cxx=%_rdict.pcm) $(TARGET_PCM) ; fi # ROOT 6

$(SUBMODULES):
	@echo "  $@ not found. Trying to initialize submodule..."
	@echo "  git submodule update --init"
	@git submodule update --init

$(BUILDDIRS):
	@exec $(MAKE) -s -C $(@:build-%=%)

$(CLEANDIRS):
	@exec $(MAKE) -s -C $(@:clean-%=%) clean

$(PREMINC): $(PREMDIR) $(PREMFILE) $(MODEL3DDIR) $(PREM3DFILE)
	@echo "#include <string>" > $@
	@echo "const std::string PREM_DIR = \"$(PREMDIR)\";" >> $@
	@echo "const std::string MODEL3D_DIR = \"$(MODEL3DDIR)\";" >> $@
	@echo "const std::string PREM_DEFAULT = \"$(PREMFILE)\";" >> $@
	@echo "const std::string PREM3D_DEFAULT = \"$(PREM3DFILE)\";" >> $@

# If running over specific models: TEST_MODELS="Model1 Model2 ..."
# Default: empty list runs all methods
TEST_MODELS ?=

# The foreach function iterates over words in TEST_MODELS,
# adding double quotes around each.
QUOTED_METHODS := $(foreach m,$(TEST_MODELS),"$(m)")

# The subst function replaces all occurrences of a space with a comma.
space := $(eval) $(eval)
comma := ,
COMMA_SEPARATED_METHODS := $(subst $(space),$(comma),$(QUOTED_METHODS))

# Enclose in curly braces
ROOT_FUNCTION_ARG := {$(COMMA_SEPARATED_METHODS)}

test: $(TARGET_LIB)
	@cd test && root -l -b -q ../tutorial/LoadOscProb.C 'TestMethods.C($(ROOT_FUNCTION_ARG))'
	@cd test && root -l -b -q ../tutorial/LoadOscProb.C 'StressTest.C($(ROOT_FUNCTION_ARG))'

clean: $(CLEANDIRS)
	@echo "  Cleaning $(PACKAGE)..."
	@rm -f $(PREMINC)
	@rm -rf tmp lib

.PHONY: $(SUBDIRS) $(BUILDDIRS) $(CLEANDIRS) test clean all
