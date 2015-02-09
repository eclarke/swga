## --- Notes ---
## Prebuilt binaries are available for Linux and Mac OS X. To use, use 
## `make prebuilt`. 
## 
## To use multithreading, pass `omp=1`. Requires gcc-4.6 or higher. Specify
## using `dsk_gcc=g++-4.6` or equivalent on Mac OS X to avoid using clang
## (default). 
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
## Using a different compiler for the dsk sources:
## $ make dsk_gcc=g++-4.9

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

# Mac OS X a) ships with an out-of-date GCC and b) symlinks `g++` to
# clang. Clang won't work with OpenMP, so to use dsk_gcc must be explicitly
# specified (e.g. dsk_gcc=g++-4.8)

uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
ifeq ($(uname_S), Darwin)  # We're on a Mac
	dsk_gcc?=/usr/bin/g++
	osx=1
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
	make -C contrib/dsk CC=$(dsk_gcc) osx=$(osx)
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


