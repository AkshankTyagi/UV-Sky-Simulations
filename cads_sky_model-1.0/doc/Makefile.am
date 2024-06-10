docfiles = skymodel_manual.tex
cadsdocdir = $(datadir)/cads/@PACKAGE@
dist_cadsdoc_DATA = $(docfiles)

if HAVE_PDFLATEX
PDFLATEX +=  --interaction=nonstopmode
dist_cadsdoc_DATA += skymodel_manual.pdf
MANNAME = skymodel_manual
MANTEXSRC = $(MANNAME).tex
MANAUX = $(MANNAME).aux
MANPDF = $(MANNAME).pdf

CLEANFILES = $(MANPDF) $(MANNAME).log $(MANNAME).idx $(MANNAME).out \
		$(MANNAME).toc $(MANAUX)

$(MANPDF): $(srcdir)/$(MANTEXSRC)
		$(PDFLATEX) $<
endif

