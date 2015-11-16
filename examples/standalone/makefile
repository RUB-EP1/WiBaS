CC = g++
CFLAGS = -Wall -ansi
SRCDIR = ../../src
SOURCES = ${SRCDIR}/phasespace_point.cc ${SRCDIR}/wibas.cc testApp.cc
OBJECTS = $(SOURCES:.cpp=.o)
INC=-I${ROOTSYS}/include -I${SRCDIR}
RLIBS = $(shell ${ROOTSYS}/bin/root-config --libs)
LDFLAGS =  ${RLIBS} -lRooFit -lRooFitCore

testApp: $(OBJECTS)
	$(CC) $(CFLAGS) $(INC) -o testApp $(OBJECTS) $(LDFLAGS)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	@echo Cleaning up...
	@rm testApp
	@echo Done.