
CXX = $(shell root-config --cxx)
CXXFLAGS = $(shell root-config --cflags) -fPIC
LDFLAGS = $(shell root-config --glibs)

VERSION = $(shell root-config --version | cut -d'.' -f1)
ifeq ($(VERSION),6)
	DICTEXE = rootcling
else
	DICTEXE = rootcint
endif

PACKAGE = MatrixDecomp
HEADERS = $(filter-out $(CURDIR)/LinkDef.h, $(wildcard $(CURDIR)/*.h))
SOURCES = $(wildcard $(CURDIR)/*.cxx)
TARGET = lib$(PACKAGE)
TARGET_LIB = $(CURDIR)/$(TARGET).so
DICTIONARY = $(CURDIR)/tmp/$(TARGET).cxx

all: $(TARGET_LIB)

$(TARGET_LIB): $(DICTIONARY) $(SOURCES)
	@echo "  Building $(PACKAGE)..."
	@$(CXX) $(CXXFLAGS) -O3 -shared -o$@ $^ $(LDFLAGS)

$(DICTIONARY): $(HEADERS) LinkDef.h
	@echo "  Making dictionary for $(PACKAGE)..."
	@mkdir -p tmp
	@$(DICTEXE) -f $@ -c -p $^
	@if [ -e $(DICTIONARY:%.cxx=%_rdict.pcm) ] ; then mv -f $(DICTIONARY:%.cxx=%_rdict.pcm) $(CURDIR)/ ; fi # ROOT 6

clean:
	@echo "  Cleaning $(PACKAGE)..."
	@rm -f $(TARGET_LIB)
	@rm -f $(TARGET_LIB:%.so=%_rdict.pcm) # ROOT 6
	@rm -f $(DICTIONARY)
	@rm -f $(DICTIONARY:%.cxx=%.h) # ROOT 5
	@rm -rf tmp

.PHONY: clean all
