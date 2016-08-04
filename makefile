CC = g++
SRCDIR = src
BINDIR = bin
INCLUDEDIR = include
TESTSRCDIR = tests
SOURCES = $(wildcard ${SRCDIR}/*cc)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cc=$(BINDIR)/%.o)
TESTSOURCES = $(wildcard ${TESTSRCDIR}/*cc)
TESTOBJECTS := $(TESTSOURCES:$(TESTSRCDIR)/%.cc=$(BINDIR)/%.o)
INC=-I${ROOTSYS}/include -I$(INCLUDEDIR)
RLIBS = $(shell ${ROOTSYS}/bin/root-config --libs)
CFLAGS = -Wall -ansi -O3 -fPIC -std=c++0x
CFLAGSEX = -Wall -ansi -O3 -std=c++0x
LDFLAGS =  ${RLIBS} -lRooFit -lRooFitCore -shared
LDFLAGSEX = ${RLIBS} -lRooFit -lRooFitCore -Wl,-rpath,./

LIBTARGET = $(BINDIR)/libwibas.so
EXAMPLEBACKGROUND = $(BINDIR)/backgroundExampleApp
EXAMPLENEREGY = $(BINDIR)/energyTestExampleApp
UNITTESTTARGET = $(BINDIR)/unitTestApp

all: $(LIBTARGET) $(EXAMPLEBACKGROUND) $(EXAMPLENEREGY) $(UNITTESTTARGET)
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

$(TESTOBJECTS) : $(BINDIR)/%.o : $(TESTSRCDIR)/%.cc
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(UNITTESTTARGET): $(TESTOBJECTS)
	$(CC) $(CFLAGSEX) $(INC) -o $@ $(TESTOBJECTS)  $(LDFLAGSEX)  -L$(BINDIR) -lwibas
	@cd $(BINDIR)/ && ../$(UNITTESTTARGET)


clean:
	@echo Cleaning up ...
	@rm -rf $(BINDIR)
	@echo Done.
