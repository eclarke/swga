VPATH = src/cliquer:swga
include src/cliquer/Makefile

SWGAHOME?=$(HOME)/.swga
SWGABIN=$(SWGAHOME)/bin
USERBASE=$(shell python -m site --user-base)

USERINST_MESSAGE=""
SWGAHOME_MESSAGE="export SWGAHOME=$(SWGAHOME)"

ifeq ($(findstr $(USERBASE)/bin, $(PATH)),)
	USERINST_MESSAGE='export PATH=$$PATH:$(USERBASE)/bin'
endif

all : cl set_finder
	pip install --editable .
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



