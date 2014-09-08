VPATH = src/cliquer:swga
include src/cliquer/Makefile

define HELP
Installation makefile for SWGA. 
Run 'make all' to compile libraries and install swga, or 'make clean' to remove any compiler intermediates.
Options:
  SWGAHOME=/some/path:       set home directory for SWGA to /some/path (default: $(SWGAHOME))
  user_install=(yes/no):     install SWGA for a single user (default: yes)
  editable_install=(yes/no): symlink repository instead of moving files (allows 'git pull' updates)
endef
export HELP

SWGAHOME?=$(HOME)/.swga
SWGABIN=$(SWGAHOME)/bin
USERBASE=$(shell python -m site --user-base)
user_install?=yes
editable_install?=yes
ifeq ($(user_install),yes)
	USER_INSTALL=--user
endif
ifeq ($(editable_install),yes)
	EDITABLE_INSTALL=--editable
endif

INSTALL_FLAGS?=$(USER_INSTALL) $(EDITABLE_INSTALL)
USERINST_MESSAGE="swga installed to $(USERBASE). Ensure $(USERBASE)/bin is in your path."
SWGAHOME_MESSAGE="Ensure \'export SWGAHOME=$(SWGAHOME)\' is in your .bashrc or .profile before running swga."
ifneq ($(USER_INSTALL),--user)
	USERINST_MESSAGE="swga installed globally"
endif

help:
	@echo "$$HELP"

all : cl set_finder
	pip install $(INSTALL_FLAGS) .
	mkdir -p $(SWGABIN)
	cp set_finder $(SWGABIN)
	cp swga/default_parameters.cfg $(SWGAHOME)/parameters.cfg
	@echo "------------------------------------"
	@echo " $(USERINST_MESSAGE)"
	@echo " $(SWGAHOME_MESSAGE)"
	@echo "------------------------------------"



