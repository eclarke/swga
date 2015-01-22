VPATH = src/cliquer:swga
include src/cliquer/Makefile

all: cl set_finder
	pip install $(opts) .
