## --- Notes ---
## Prebuilt binaries are available for Linux and Mac OS X. To compile from 
## source instead, use `make source`.
## 
## For Mac OS X, pass `osx=1`.
## Yosemite, gcc 4.8 or higher may need to be installed.
##
## For a user-only install (no root needed), pass `user=1`.
## This is incompatible with virtualenvs or with homebrew python on Mac OS X.
##
## For an editable install (devs only!) pass `editable=1`
## This allows changes made to this directory to change SWGA's behavior.
##
## --- Examples ---
## Installing on Linux using prebuilt binaries:
## $ make user=1
##
## Installing on Mac OS X from source:
## $ make compile user=1 osx=1
##
## Using a different compiler for the dsk sources:
## $ make compile user=1 osx=1 dsk_gcc=g++-4.9

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
	dsk_gcc?=g++-4.8
	binaries=contrib/bin/osx
else
	dsk_gcc?=g++
	binaries=contrib/bin/linux
endif


all: compile


pip:
	pip --version > /dev/null

compile: pip
	make -C contrib/cliquer
	make -C contrib/dsk CC=$(dsk_gcc)
	mkdir -p swga/bin
	cp contrib/cliquer/set_finder swga/bin/
	cp contrib/dsk/dsk swga/bin/
	cp contrib/dsk/parse_results swga/bin/
	pip install $(pip_flags) .
ifneq ($(swga_cmd_path),)
	python swga/data/finished_message.py "$(swga_cmd_path)"
endif


prebuilt: pip
	mkdir -p swga/bin
	cp $(binaries)/* swga/bin/
	pip install $(pip_flags) .
ifneq ($(swga_cmd_path),)
	python swga/data/finished_message.py "$(swga_cmd_path)"
endif


