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

USERINST_MESSAGE=""
SWGAHOME_MESSAGE="export SWGAHOME=$(SWGAHOME)"

ifeq ($(findstr $(USERBASE)/bin, $(PATH)),)
	USERINST_MESSAGE='export PATH=$$PATH:$(USERBASE)/bin'
endif

help:
	@echo "$$HELP"

all : cl set_finder
	pip install --user --editable .
	mkdir -p $(SWGABIN)
	cp set_finder $(SWGABIN)
	cp swga/default_parameters.cfg $(SWGAHOME)/default_parameters.cfg
	@echo "------------------------------------"
	@echo "Install succeeded! "
	@echo "Before using SWGA, run the following commands or put them in your \n\
	shell config file (such as .bashrc, .profile, or .bash_profile). \n\
	If not placed in your shell config file, these commands must be run each \n\
	time before using SWGA:"
	@echo
	@echo   $(SWGAHOME_MESSAGE)
	@echo   $(USERINST_MESSAGE)
	@echo "------------------------------------"



