
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

PACKAGE = MatrixDecomp
SOURCES := $(wildcard $(CURDIR)/*.cxx)
TARGET = lib$(PACKAGE)
TARGET_LIB = $(CURDIR)/$(TARGET).so

all: $(TARGET_LIB)

$(TARGET_LIB): $(SOURCES)
	@echo "  Building $(PACKAGE)..."
	@$(CXX) $(CXXFLAGS) -O3 -shared -o$@ $^ $(LDFLAGS)

clean:
	@echo "  Cleaning $(PACKAGE)..."
	@rm -f $(TARGET_LIB)

.PHONY: clean all
