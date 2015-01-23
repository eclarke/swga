## --- Notes ---
## For Mac OS X, pass `osx=1` and `CC=g++-4.8` (or higher).
## gcc 4.8 or higher may need to be installed.
##
## For a user-only install (no root needed), pass `user=1`.
## This is incompatible with virtualenvs or with homebrew python on Mac OS X.
##
## For an editable install (devs only!) pass `editable=1`
## This allows changes made to this directory to change SWGA's behavior.
##
## --- Example ---
## A Mac OS X user install:
## make user=1 osx=1 CC=g++-4.8

ifeq ($(user),1)
	pip_flags+= --user
endif
user_base = $(shell python -m site --user-base)	

ifeq ($(findstring $(user_base)/bin,$(PATH)),)
	swga_cmd_path=$(user_base)/bin
endif


ifeq ($(editable),1)
	pip_flags+= --editable
endif

ifeq ($(osx),1)
	binaries=contrib/bin/osx
else
	binaries=contrib/bin/linux
endif

## DEBUG
swga_cmd_path = $(user_base)/blah/

all:
	make -C contrib/cliquer
	make -C contrib/dsk
	mkdir -p swga/bin
	cp contrib/cliquer/set_finder swga/bin/
	cp contrib/dsk/dsk swga/bin/
	cp contrib/dsk/parse_results swga/bin/
	pip install $(pip_flags) .

ifneq ($(swga_cmd_path),)
	python swga/data/finished_message.py "$(swga_cmd_path)"
endif

prebuilt:
	mkdir -p swga/bin
	cp $(binaries)/* swga/bin/
	pip install $(pip_flags) .
ifeq ($(swga_cmd_path),)
	python swga/data/finished_message.py "$(swga_cmd_path)"
endif
