bin_PROGRAMS = sky_model
sky_model_SOURCES = sky_model.c calculations.c coordinates.c gal_bkgd.c \
                         parse_par.c read_write.c stars.c zodiacal_programs.c

noinst_HEADERS = parse_par.h sky_model.h
cadsdatadir = $(datadir)/cads/data
AM_CFLAGS = -DCADSDATA=\"$(cadsdatadir)/\" -lm