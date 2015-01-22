VPATH = src/cliquer:swga
include src/cliquer/Makefile

in_venv = $(shell python -c "import sys; print(hasattr(sys, 'real_prefix'))")
user_base = $(shell python -m site --user-base)
cmd_prefix =

ifeq ($(in_venv), False)
	pip_opts = --user
	ifeq ($(findstring $(user_base)/bin,$(PATH)),)
		cmd_prefix = $(user_base)/bin/
	endif
endif

opts?=
pip_opts += $(opts)


all: cl set_finder
	@echo $(PATH)
	pip install $(pip_opts) .
	python swga/data/finished_message.py "$(cmd_prefix)"
