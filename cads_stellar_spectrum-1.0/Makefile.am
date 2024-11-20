# not a GNU package. You can remove this line, if
# have all needed files, that a GNU package needs
AUTOMAKE_OPTIONS = foreign 1.4

SUBDIRS = src
cadsdocdir = $(datadir)/cads/@PACKAGE@
dist_cadsdoc_DATA = README.txt AUTHORS.txt COPYING.txt ChangeLog.txt INSTALL.txt

# if you change cadsdatadir, update it in src/Makefile.am as well
# -- I've no idea about passing VAR from one Makefile.am to another :(
cadsdatadir = $(datadir)/cads/data

# incase the following .mk files need to be generated again,
# try the following from this directory level:
#-- echo -n "KURUCZ_FILES =" > spectral_data/kurucz_files.mk
#-- find spectral_data/KURUCZ_FLUX -type f -print | \
#-- sed 's/^/  /;$q;s/$/ \\/' >> spectral_data/kurucz_files.mk
#--
#-- echo -n "SPDATA_FILES =" > spectral_data/spdata_files.mk
#-- find spectral_data/MKSpDATA -type f -print | \
#-- sed 's/^/  /;$q;s/$/ \\/' >> spectral_data/spdata_files.mk
#
# thanks to a question by Jeremy Salwen and adl's response on stackoverflow
# --Reks
include spectral_data/kurucz_files.mk
include spectral_data/spdata_files.mk

spdatadir = $(cadsdatadir)/spectral_data
dist_spdata_DATA = spectral_data/crossec0.dat spectral_data/wave.dat

kuruczdir = $(spdatadir)/KURUCZ_FLUX/P00
dist_kurucz_DATA = $(KURUCZ_FILES)

mkspdir = $(spdatadir)/MKSpDATA
dist_mksp_DATA = $(SPDATA_FILES)

