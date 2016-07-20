CC = g++
SRCDIR = src
BINDIR = bin
INCLUDEDIR = include
SOURCES = $(wildcard ${SRCDIR}/*cc)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cc=$(BINDIR)/%.o)
INC=-I${ROOTSYS}/include -I$(INCLUDEDIR)
RLIBS = $(shell ${ROOTSYS}/bin/root-config --libs)
CFLAGS = -Wall -ansi -O3 -fPIC
CFLAGSEX = -Wall -ansi -O3
LDFLAGS =  ${RLIBS} -lRooFit -lRooFitCore -shared
LDFLAGSEX = ${RLIBS} -lRooFit -lRooFitCore -Wl,-rpath,./

LIBTARGET = $(BINDIR)/libwibas.so
EXAMPLETARGET = $(BINDIR)/testApp

all: $(LIBTARGET) $(EXAMPLETARGET)
	@mkdir -p bin

$(OBJECTS): $(BINDIR)/%.o : $(SRCDIR)/%.cc
	@mkdir -p $(@D) 
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(LIBTARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

$(EXAMPLETARGET):  examples/standalone/testApp.cc
	$(CC)  $(CFLAGSEX) $(INC)  -o $@ $< $(LDFLAGSEX)  -L$(BINDIR) -lwibas

clean:
	@echo Cleaning up ...
	@rm -rf $(BINDIR)
	@echo Done.
