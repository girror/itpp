dnl @synopsis ACX_CBLAS
dnl @author Adam Piatyszek <ediap@users.sourceforge.net>
dnl @version 2006-01-31
dnl
dnl This macro looks for a library that implements the CBLAS
dnl interface for BLAS (see http://www.netlib.org/blas/). On
dnl success, it sets the CBLAS_LIBS output variable to hold the
dnl requisite library linkages. Besides, it defines HAVE_CBLAS.
dnl
dnl To link with CBLAS, you should link with:
dnl
dnl     $CBLAS_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically. FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS), and
dnl is sometimes necessary in order to link with F77 libraries. Users
dnl will also need to use AC_F77_DUMMY_MAIN (see the autoconf manual),
dnl for the same reason.
dnl
dnl The user may also use --with-cblas=<lib> in order to use some
dnl specific CBLAS library <lib>. In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as was
dnl used to compile the BLAS library.
dnl
dnl @license GPLWithACException

AC_DEFUN([ACX_CBLAS], [
AC_REQUIRE([ACX_BLAS])

AC_ARG_WITH(cblas,
  [AC_HELP_STRING([--with-cblas=<lib>], [use CBLAS library <lib>])])
case $with_cblas in
  yes | "") ;;
  no) acx_cblas_ok=disabled ;;
  -* | */* | *.a | *.so | *.so.* | *.o) CBLAS_LIBS="$with_cblas" ;;
  *) CBLAS_LIBS="-l$with_cblas" ;;
esac

test x"$acx_cblas_ok" != xdisabled && acx_cblas_ok=no

# We cannot use CBLAS if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
  acx_cblas_ok=noblas
fi

# First, check CBLAS_LIBS environment variable
if test "x$CBLAS_LIBS" != x; then
  save_LIBS="$LIBS"; LIBS="$CBLAS_LIBS $BLAS_LIBS $LIBS $FLIBS"
  AC_MSG_CHECKING([for cblas_sgemm in $CBLAS_LIBS])
  AC_TRY_LINK_FUNC(cblas_sgemm, [acx_cblas_ok=yes], [CBLAS_LIBS=""])
  AC_MSG_RESULT($acx_cblas_ok)
  LIBS="$save_LIBS"
fi

# CBLAS linked to by default?  (is sometimes included in BLAS lib)
if test $acx_cblas_ok = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS"
  AC_CHECK_FUNC(cblas_sgemm, [acx_cblas_ok=yes])
  LIBS="$save_LIBS"
fi

# CBLAS linked to by default using FLIBS?
if test $acx_cblas_ok = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
  AC_CHECK_FUNC(cblas_sgemm, [acx_cblas_ok=yes])
  LIBS="$save_LIBS"
fi

# CBLAS from ACML linked to by default using FLIBS?
cblas_acml_ok=no
if test $acx_cblas_ok = no; then
  save_LIBS="$LIBS"; LIBS="$LIBS $BLAS_LIBS $FLIBS"
  AC_CHECK_FUNC(sgemm, [acx_cblas_ok=yes; cblas_acml_ok=yes])
  LIBS="$save_LIBS"
fi

# Generic CBLAS library?
for cblas in cblas gslcblas; do
  if test $acx_cblas_ok = no; then
    save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS $FLIBS"
    AC_CHECK_LIB($cblas, cblas_sgemm, [acx_cblas_ok=yes; CBLAS_LIBS="-l$cblas"])
    LIBS="$save_LIBS"
  fi
done

AC_SUBST(CBLAS_LIBS)

# Finally, define HAVE_CBLAS
if test x"$acx_cblas_ok" = xyes; then
  AC_DEFINE(HAVE_CBLAS, 1, [Define if you have CBLAS library.])
  if test x"$cblas_acml_ok" = xyes; then
    AC_DEFINE(HAVE_CBLAS_ACML, 1, [Define if you have ACML CBLAS library.])
  fi
fi

])dnl ACX_CBLAS
