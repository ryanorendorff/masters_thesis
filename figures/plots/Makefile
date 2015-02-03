# This is to make the stupid expansion in clean work.
SHELL=/bin/bash

SUBDIRS:=$(shell ls -d */)

#######################################################################
#                               Targets                               #
#######################################################################
.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@
