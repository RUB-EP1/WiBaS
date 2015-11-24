CC = g++
CFLAGS = -Wall -ansi -O3
SRCDIR = ../../src
SOURCES = $(wildcard ${SRCDIR}/*cc) testApp.cc
OBJECTS = $(SOURCES:.cc=.o)
INC=-I${ROOTSYS}/include -I${SRCDIR}
RLIBS = $(shell ${ROOTSYS}/bin/root-config --libs)
LDFLAGS =  ${RLIBS} -lRooFit -lRooFitCore

%.o: %.cc
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

testApp: $(OBJECTS)
	$(CC) $(CFLAGS) $(INC) -o testApp $(OBJECTS) $(LDFLAGS)

clean:
	@echo Cleaning up...
	@rm *.o
	@rm ${SRCDIR}/*o
	@rm testApp
	@echo Done.