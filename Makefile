VPATH = src/cliquer:swga
include src/cliquer/Makefile
SWGAHOME?=$(HOME)/.swga
SWGABIN=$(SWGAHOME)/bin
USERBASE=$(shell python -m site --user-base)
USER_INSTALL?=--user
EDITABLE_INSTALL?=--editable
INSTALL_FLAGS?=$(USER_INSTALL) $(EDITABLE_INSTALL)
USERINST_MESSAGE="swga installed to $(USERBASE). Ensure $(USERBASE)/bin is in your path.\n"
SWGAHOME_MESSAGE="Ensure 'export SWGAHOME=$(SWGAHOME)' is in your .bashrc or .profile before running swga."
ifneq ($(USER_INSTALL),--user)
	MESSAGE="swga installed globally"
endif
all : cl set_finder
	pip install $(INSTALL_FLAGS) .
	mkdir -p $(SWGABIN)
	cp set_finder $(SWGABIN)
	cp swga/default_parameters.cfg $(SWGAHOME)/parameters.cfg
	@echo "------------------------------------"
	@echo " $(USERINST_MESSAGE)"
	@echo " $(SWGAHOME_MESSAGE)"
	@echo "------------------------------------"



