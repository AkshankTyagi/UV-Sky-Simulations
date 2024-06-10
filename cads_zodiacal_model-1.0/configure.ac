AC_INIT([CADS_zodiacal_model], [1.0], [cads@iiap.res.in])
AC_LANG_C
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdio.h stdlib.h string.h math.h])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC

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


AC_OUTPUT(Makefile src/Makefile)

echo
echo " +------------------------------------------------------------------+"
echo " | Configuration complete                                           |"
echo " | Type 'make' to compile the software, followed by 'make install'  |"
echo " |                   All the very best! :-)                         |"
echo " |                                             --CADS Software Team |"
echo " +------------------------------------------------------------------+"
echo

