CC = g++
SRCDIR = src
BINDIR = bin
INCLUDEDIR = include
SOURCES = $(wildcard ${SRCDIR}/*cc)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cc=$(BINDIR)/%.o)
INC=-I${ROOTSYS}/include -I$(INCLUDEDIR)
RLIBS = $(shell ${ROOTSYS}/bin/root-config --libs)
CFLAGS = -Wall -ansi -O3 -fPIC -std=c++0x
CFLAGSEX = -Wall -ansi -O3 -std=c++0x
LDFLAGS =  ${RLIBS} -lRooFit -lRooFitCore -shared
LDFLAGSEX = ${RLIBS} -lRooFit -lRooFitCore -Wl,-rpath,./

LIBTARGET = $(BINDIR)/libwibas.so
EXAMPLEBACKGROUND = $(BINDIR)/backgroundExampleApp
EXAMPLENEREGY = $(BINDIR)/energyTestExampleApp

all: $(LIBTARGET) $(EXAMPLEBACKGROUND) $(EXAMPLENEREGY)
	@mkdir -p bin

$(OBJECTS): $(BINDIR)/%.o : $(SRCDIR)/%.cc
	@mkdir -p $(@D) 
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(LIBTARGET) : $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

$(EXAMPLEBACKGROUND):  examples/standalone/backgroundExampleApp.cc
	$(CC)  $(CFLAGSEX) $(INC)  -o $@ $< $(LDFLAGSEX)  -L$(BINDIR) -lwibas

$(EXAMPLENEREGY):  examples/standalone/energyTestExampleApp.cc
	$(CC) $(CFLAGSEX) $(INC) -o $@ $< $(LDFLAGSEX)  -L$(BINDIR) -lwibas

clean:
	@echo Cleaning up ...
	@rm -rf $(BINDIR)
	@echo Done.
