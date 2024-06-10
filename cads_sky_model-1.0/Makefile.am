# not a GNU package. You can remove this line, if
# have all needed files, that a GNU package needs
AUTOMAKE_OPTIONS = foreign 1.4

SUBDIRS = src doc
cadsdocdir = $(datadir)/cads/@PACKAGE@
dist_cadsdoc_DATA = README.txt AUTHORS.txt COPYING.txt ChangeLog.txt INSTALL.txt

# if you change cadsdatadir, update it in src/Makefile.am as well
# -- I've no idea about passing VAR from one Makefile.am to another :(
cadsdatadir = $(datadir)/cads/data

stdatadir = $(cadsdatadir)/stellar_data
dfdatadir = $(cadsdatadir)/diffuse_data
zodatadir = $(cadsdatadir)/zodiacal_data
fldatadir = $(cadsdatadir)/sample_filters

dist_stdata_DATA = data/stellar_data/cross_sec.txt      \
                   data/stellar_data/hip_main.txt
dist_dfdata_DATA = data/diffuse_data/bkgd_spectrum.txt  \
                   data/diffuse_data/fuv_allsky.txt     \
                   data/diffuse_data/nuv_allsky.txt
dist_zodata_DATA = data/zodiacal_data/zodiacal_spec.txt \
                   data/zodiacal_data/leinert_dist.txt
dist_fldata_DATA = data/sample_filters/galex_FUV.txt    \
                   data/sample_filters/galex_NUV.txt    \
                   data/sample_filters/uvit_FUVBAF2.txt \
                   data/sample_filters/uvit_NUVN2.txt

include data/stellar_data/castelli_files.mk
castellidir = $(stdatadir)/castelli
dist_castelli_DATA = $(CASTELLI_FILES)

