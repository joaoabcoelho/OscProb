CXX = $(shell root-config --cxx)
CXXFLAGS = $(shell root-config --cflags) -fPIC
LDFLAGS = $(shell root-config --glibs)

PREFIX = ../..
BINDIR = bin

override LDFLAGS += -L$(PREFIX)/lib $(PREFIX)/lib/libOscProb.so

INCDIRS = -I$(PREFIX)/inc -Iinc -I$(PREFIX)/eigen

SOURCE = $(wildcard src/*.cc)
TARGET = $(patsubst %.cc,%,$(shell basename $(SOURCE)))

$(BINDIR)/$(TARGET): $(SOURCE) $(PREFIX)/lib/libOscProb.so
	@mkdir -p ${BINDIR}
	@echo "Compiling ${TARGET} from ${SOURCE}..."
	@$(CXX) -O3 -o $@ $< $(LDFLAGS) $(CXXFLAGS) $(INCDIRS)

phony: $(PREFIX)/lib/libOscProb.so

clean:
	@rm -f $(BINDIR)/*
