AC_INIT([CADS_sky_model], [1.0], [cads@iiap.res.in])

AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_STDC
AC_LANG_C
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
CFLAGS="-g"

#-- check if we have pdf-latex to make user manual.
AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
if test -z "$PDFLATEX"; then
AC_MSG_WARN([Wont make PDF version of the user manual])
fi
AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")
#---------------------------------------------------------------------#
#           Now we check for cfitsio headers and library              #
#---------------------------------------------------------------------#
LIBS="$LIBS -L/usr/lib64 -L/usr/lib -L/usr/local/lib -L/usr/local/lib64\
 -lm -lcfitsio"
CFLAGS="$CFLAGS -I/usr/include -I/usr/local/include\
 -I/usr/include/cfitsio -I/usr/local/include/cfitsio"

fitserror="
 +------------------------------------------------------------------+
 | Unable to find cfitsio distribution. You may need to get it from |
 |   http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html  |
 +------------------------------------------------------------------+
"
AC_ARG_WITH(cfitsio,
     [ --with-cfitsio[=DIR]     Location of the cfitsio distribution],
     [CFITS=$with_cfitsio]
)

if test x"$CFITS" != "x"; then
    LIBS="-L$CFITS/lib64 -L$CFITS/lib -L$CFITS $LIBS"
    CFLAGS="$CFLAGS -I$CFITS -I$CFITS/include -I$CFITS/include/cfitsio"
fi

AC_MSG_CHECKING(for cfitsio)
AC_LINK_IFELSE([AC_LANG_PROGRAM([[
    #include <stdio.h>
    #include <stdlib.h>
    #include <fitsio2.h>]],
    [[    fitsfile *fptr]])],
    [cfitsio=yes],[cfitsio=no])
AC_MSG_RESULT([$cfitsio])

if test "$cfitsio" != "yes"; then
    AC_MSG_ERROR($fitserror)
fi

#---------------------------------------------------------------------#
#           Now we check for libnova headers and library              #
#---------------------------------------------------------------------#
LIBS="$LIBS -L/usr/lib -L/usr/lib64 -L/usr/local/lib \
-L/usr/local/lib64 -lnova"
CFLAGS="$CFLAGS -I/usr/include -I/usr/include/libnova \
-I/usr/local/include -I/usr/local/include/libnova"

lnovaerror="
 +------------------------------------------------------------------+
 | Unable to find libnova distribution. You may need to get it from |
 |                  http://libnova.sourceforge.net/                 |
 +------------------------------------------------------------------+
"
AC_ARG_WITH(libnova,
     [ --with-libnova[=DIR]     Location of the libnova distribution],
     [LNOVA=$with_libnova]
)

if test x"$LNOVA" != 'x'; then
    LIBS="-L$LNOVA/lib64 -L$LNOVA/lib -L$LNOVA $LIBS"
    CFLAGS="$CFLAGS -I$LNOVA -I$LNOVA/include -I$LNOVA/include/libnova"
fi

AC_MSG_CHECKING(for libnova)
AC_LINK_IFELSE([AC_LANG_PROGRAM([[
    #include <stdio.h>
    #include <stdlib.h>
    #include <libnova/libnova.h>]],
    [[    struct ln_lnlat_posn observer;]])],
    [libnova=yes],[libnova=no])
AC_MSG_RESULT([$libnova])

if test "$libnova" != "yes"; then
    AC_MSG_ERROR($lnovaerror)
fi

AC_OUTPUT(Makefile src/Makefile doc/Makefile)
echo
echo " +------------------------------------------------------------------+"
echo " | Configuration complete                                           |"
echo " | Type 'make' to compile the software, followed by 'make install'  |"
echo " |                   All the very best! :-)                         |"
echo " |                                             --CADS Software Team |"
echo " +------------------------------------------------------------------+"
echo
