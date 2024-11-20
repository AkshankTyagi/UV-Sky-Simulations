AC_INIT([CADS_stellar_spectrum], [1.0], [cads@iiap.res.in])
AC_LANG_C
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for header files.
AC_CHECK_HEADERS([stdio.h stdlib.h string.h math.h])

# Checks for library functions.
AC_FUNC_MALLOC

AC_OUTPUT(Makefile src/Makefile)

echo
echo " +------------------------------------------------------------------+"
echo " | Configuration complete                                           |"
echo " | Type 'make' to compile the software, followed by 'make install'  |"
echo " |                   All the very best! :-)                         |"
echo " |                                             --CADS Software Team |"
echo " +------------------------------------------------------------------+"
echo
