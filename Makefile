# This is to make the stupid expansion in clean work.
SHELL=/bin/bash

SUBDIRS=figures/plots

DOCUMENT=mphil_cam_orendorff

TEX=pdflatex
BIBTEX=bibtex

#######################################################################
#                               Targets                               #
#######################################################################
.PHONY: subdirs $(SUBDIRS)

${DOCUMENT}.pdf: subdirs ${DOCUMENT}.tex
	${TEX} ${DOCUMENT}
	${BIBTEX} ${DOCUMENT}
	${TEX} ${DOCUMENT}
	${TEX} ${DOCUMENT}

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@


