bin_PROGRAMS = uvs_stellar_spectrum

# the library search path.
uvs_stellar_spectrum_SOURCES = signal_handler.c stellar_spectrum.c getopt.c kurucz.c \
read_write.c
noinst_HEADERS = stellar_spectrum.h
cadsdatadir = $(datadir)/cads/data

AM_CFLAGS = -DCADSDATA=\"$(cadsdatadir)/\" -lm