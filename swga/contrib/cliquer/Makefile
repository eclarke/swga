
##### Configurable options:

## Compiler:
CC=gcc
#CC=cc

## Compiler flags:

# GCC:  (also -march=pentium etc, for machine-dependent optimizing)
CFLAGS=-Wall -O3 -fomit-frame-pointer -funroll-loops

# GCC w/ debugging:
# CFLAGS=-Wall -g -DINLINE=

# Compaq C / Digital C:
#CFLAGS=-arch=host -fast

# SunOS:
#CFLAGS=-fast

## Program options:

# Enable long options for cl (eg. "cl --help"), comment out to disable.
# Requires getopt_long(3)  (a GNU extension)
LONGOPTS = -DENABLE_LONG_OPTIONS


##### End of configurable options


all: cl set_finder


testcases: testcases.o cliquer.o graph.o reorder.o
	$(CC) $(LDFLAGS) -o $@ testcases.o cliquer.o graph.o reorder.o

cl: cl.o cliquer.o graph.o reorder.o
	$(CC) $(LDFLAGS) -o $@ cl.o cliquer.o graph.o reorder.o

set_finder: set_finder.o cliquer.o graph.o reorder.o
	$(CC) $(LDFLAGS) -o $@ set_finder.o cliquer.o graph.o reorder.o

cl.o testcases.o cliquer.o graph.o reorder.o: cliquer.h set.h graph.h misc.h reorder.h Makefile cliquerconf.h

cl.o: cl.c
	$(CC) $(CFLAGS) $(LONGOPTS) -o $@ -c $<

set_finder.o: set_finder.c
	$(CC) $(CFLAGS) $(LONGOPTS) -o $@ -c $<

clean:
	rm -f *.o *~ cl testcases set_finder

backup:
	mkdir "`date "+backup-%Y-%m-%d-%H-%M"`" 2>/dev/null || true
	cp * "`date "+backup-%Y-%m-%d-%H-%M"`"  2>/dev/null || true

test: testcases
	./testcases
