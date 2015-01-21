VPATH = src/cliquer:swga
include src/cliquer/Makefile

USERINST_MESSAGE=\nBefore using SWGA, you need to run the file activate_swga.sh. \nThis file needs to be run every time you start a new session or terminal window before using SWGA. \nAlternatively, you can copy its contents into your ~/.bashrc or ~/.profile.

# Tests whether or not we're installing into a virtualenv- if so,
# we omit the --user part of the pip install command (since they're incompatible)
PIP_VERSION=$(shell pip -V | cut -f2 -d" ")
IN_VIRTUALENV=$(shell python -c "import sys; print(hasattr(sys, 'real_prefix'))")
PIP_OPTIONS=--editable
# If the user is root, then we do not install an editable copy (nor do we
# install it to the user site-packages)
ifeq ($(USER), root)
	PIP_OPTIONS=
else
	ifeq ($(IN_VIRTUALENV), False)
		PIP_OPTIONS=--user --editable
		USERBASE=$(shell python -m site --user-base)
		ifeq ($(findstr $(USERBASE)/bin, $(PATH)),)
			USERINST_SCRIPT='export PATH=$$PATH:$(USERBASE)/bin'
		endif
	endif
endif

all : cl set_finder
	pip install $(PIP_OPTIONS) .
	echo $(USERINST_SCRIPT) > activate_swga.sh
	chmod +x activate_swga.sh
	@echo "------------------------------------"
	@echo "------------------------------------"
	@echo ""
	@echo "Install succeeded! "
	@echo "$(USERINST_MESSAGE)"
	@echo ""
	@echo "In a different directory, run 'init_swga' to set up a new workspace and get started."
	@echo ""
	@echo "------------------------------------"
	@echo "------------------------------------"



