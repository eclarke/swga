VPATH = src/cliquer:swga
include src/cliquer/Makefile

in_venv = $(shell python -c "import sys; print(hasattr(sys, 'real_prefix'))")
user_base = $(shell python -m site --user-base)
cmd_prefix =
ifneq ($(findstring $(opts),--user),)
	ifeq ($(findstring $(user_base)/bin,$(PATH)),)
		cmd_prefix = $(user_base)/bin/
	endif
endif

all: cl set_finder
	pip install $(opts) .
	python swga/data/finished_message.py "$(cmd_prefix)"
