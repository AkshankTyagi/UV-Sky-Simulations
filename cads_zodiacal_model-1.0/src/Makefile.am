bin_PROGRAMS = uvs_zodiacal_model

# the library search path.
uvs_zodiacal_model_SOURCES = signal_handler.c zodiacal_model.c read_write.c
noinst_HEADERS = zodiacal_model.h

cadsdatadir = $(datadir)/cads/data
dist_cadsdata_DATA = leinert_dist.txt zodiacal_spec.txt
AM_CFLAGS = -DCADSDATA=\"$(cadsdatadir)/\"