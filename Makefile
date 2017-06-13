
CXX = g++
CXXFLAGS = $(shell root-config --cflags) -fPIC
LDFLAGS = $(shell root-config --glibs)

PACKAGE = OscProb
HEADERS = $(filter-out $(CURDIR)/LinkDef.h, $(wildcard $(CURDIR)/*.h))
SOURCES = $(wildcard $(CURDIR)/*.cxx)
TARGET_LIB = $(CURDIR)/lib$(PACKAGE).so
DICTIONARY = $(CURDIR)/tmp/dic$(PACKAGE).cxx

# GSL library
GSL_INCS = $(shell gsl-config --cflags)
GSL_LIBS = $(shell gsl-config --libs)

override CXXFLAGS += $(GSL_INCS)
override LDFLAGS  += $(GSL_LIBS)

# the sets of directories to do various things in
MATRIX=$(CURDIR)/MatrixDecomp/libMatrixDecomp.so
SUBDIRS = MatrixDecomp
BUILDDIRS = $(SUBDIRS:%=build-%)
CLEANDIRS = $(SUBDIRS:%=clean-%)

PREMDIR = $(CURDIR)/PremTables
PREMFILE = $(PREMDIR)/prem_default.txt
PREMINC = $(CURDIR)/prem_default.hpp

all: $(BUILDDIRS) $(PREMINC) $(TARGET_LIB)

$(TARGET_LIB): $(DICTIONARY) $(SOURCES) $(MATRIX)
	@echo "  Building $(PACKAGE)..."
	@g++ -shared -O3 -o$@ $(LDFLAGS) $(CXXFLAGS) -I$(ROOTSYS)/include $^

$(DICTIONARY): $(HEADERS) LinkDef.h
	@echo "  Making dictionary for $(PACKAGE)..."
	@mkdir -p tmp
	@rootcint -f $@ -c $(GSL_INCS) -p $^

$(BUILDDIRS):
	@exec $(MAKE) -s -C $(@:build-%=%)

$(CLEANDIRS):
	@exec $(MAKE) -s -C $(@:clean-%=%) clean

$(PREMINC): $(PREMDIR) $(PREMFILE)
	@echo "#include <string>" > $@
	@echo "const std::string PREM_DIR = \"$(PREMDIR)\";" >> $@
	@echo "const std::string PREM_DEFAULT = \"$(PREMFILE)\";" >> $@

clean: $(CLEANDIRS)
	@echo "  Cleaning $(PACKAGE)..."
	@rm -f $(TARGET_LIB)
	@rm -f $(PREMINC)
	@rm -f $(DICTIONARY)
	@rm -f $(DICTIONARY:%.cxx=%.h)
	@rm -rf tmp


.PHONY: $(SUBDIRS) $(BUILDDIRS) $(CLEANDIRS) clean all

