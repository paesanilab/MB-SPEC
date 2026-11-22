dnl (c) countless GNU projects  

AC_PREREQ(2.54)

dnl define the version of this package
m4_define([version_major], [0])
m4_define([version_minor], [1])
m4_define([version_micro], [0])

m4_define([the_version],[version_major.version_minor.version_micro])

dnl initialize autoconf
AC_INIT([mb-spec], [the_version])
AM_INIT_AUTOMAKE([foreign 1.9 no-define dist-bzip2 no-dist-gzip])

AC_CONFIG_SRCDIR([src/tools/io-tools.cpp])
AC_CONFIG_HEADERS(config.h)

dnl support silent build rules
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

dnl need this one
AC_PROG_CC
AC_PROG_CXX
AC_LANG([C++])

AM_MAINTAINER_MODE

AC_ARG_ENABLE(ansi,
    AC_HELP_STRING([--enable-ansi=@<:@no/yes@:>@],
                   [pass -ansi -pedantic to gcc @<:@default=no@:>@]),,
                   enable_ansi=no)

dnl build time sanity check
AM_SANITY_CHECK

dnl checks for programs
AC_PROG_RANLIB
AC_PROG_INSTALL
AC_PROG_MAKE_SET

changequote(,)dnl
if test "x$GXX" = "xyes"; then
    case " $CXXFLAGS " in
        *[\ \	]-Wall[\ \	]*) ;;
        *) CXXFLAGS="$CXXFLAGS -Wall" ;;
    esac

    if test "x$enable_ansi" = "xyes"; then
        case " $CXXFLAGS " in
            *[\ \	]-ansi[\ \	]*) ;;
            *) CXXFLAGS="$CXXFLAGS -ansi" ;;
        esac

        case " $CXXFLAGS " in
            *[\ \	]-pedantic[\ \	]*) ;;
            *) CXXFLAGS="$CXXFLAGS -pedantic" ;;
        esac
    fi
fi
changequote([,])dnl

dnl GSL
saved_LIBS="$LIBS"
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],[],
[AC_MSG_ERROR([
************************************************************************
**  cannot continue without GSL (http://www.google.com/search?q=GSL)  **
************************************************************************])])
LIBS="$saved_LIBS"

GSL_LIBS="-lgsl -lgslcblas"
AC_SUBST(GSL_LIBS)

dnl MKL is required by the matrix utilities.
AC_ARG_ENABLE([mkl],
    AC_HELP_STRING([--enable-mkl=@<:@yes/no@:>@],
                   [use required MKL LAPACKE routines @<:@default=yes; no is unsupported@:>@]),,
                   enable_mkl=yes)

if test "x$enable_mkl" != "xyes"; then
    AC_MSG_ERROR([MKL is required by MB-SPEC; --disable-mkl is unsupported.
Install Intel MKL and set --with-mkl-cflags and --with-mkl-libs for your installation.])
fi

AC_ARG_WITH([mkl-cflags],
    AC_HELP_STRING([--with-mkl-cflags=FLAGS],
                   [MKL compile flags]),
                   [MKL_CFLAGS=$with_mkl_cflags],
                   [MKL_CFLAGS='-m64 -I/opt/intel/mkl/include'])

AC_ARG_WITH([mkl-libs],
    AC_HELP_STRING([--with-mkl-libs=FLAGS],
                   [MKL linking flags]),
                   [MKL_LIBS=$with_mkl_libs],
                   [MKL_LIBS='-L/opt/intel/lib/intel64 -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread'])

saved_LIBS=$LIBS
saved_CXXFLAGS=$CXXFLAGS
LIBS="$LIBS $MKL_LIBS"
CXXFLAGS="$CXXFLAGS $MKL_CFLAGS"
AC_MSG_CHECKING([for required MKL LAPACKE routines])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include <mkl_lapacke.h>]], [[
        double a[1] = {1.0};
        double eigenvalues[1];
        lapack_int pivots[1];
        lapack_int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 1, 1, a, 1, pivots);
        info += LAPACKE_dgetri(LAPACK_ROW_MAJOR, 1, a, 1, pivots);
        info += LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', 1, a, 1, eigenvalues);
        return info;
    ]])],
    [mkl_usable=yes],
    [mkl_usable=no])
CXXFLAGS=$saved_CXXFLAGS
LIBS=$saved_LIBS
AC_MSG_RESULT([$mkl_usable])
if test "x$mkl_usable" != "xyes"; then
    AC_MSG_ERROR([MKL is required, but mkl_lapacke.h or the LAPACKE_dgetrf,
LAPACKE_dgetri, and LAPACKE_dsyev routines could not be compiled and linked.
Set --with-mkl-cflags and --with-mkl-libs for your Intel MKL installation.
See config.log for the compiler or linker error.])
fi

AC_DEFINE([USE_MKL], [1], [use required MKL LAPACKE routines])
AC_SUBST(MKL_CFLAGS)
AC_SUBST(MKL_LIBS)

dnl fftw3 needed for pimd
dnl PKG_CHECK_MODULES([FFTW], [fftw3])
dnl

AC_CONFIG_FILES([
  Makefile
  src/Makefile
  src/tools/Makefile
  src/infrared/Makefile
  src/raman/Makefile
  src/sfg/Makefile
  src/power/Makefile
  src/unittests/Makefile
])

AC_OUTPUT

dnl
dnl summary
dnl

echo "<=============================================================>
 Configuration summary :
<------------------------------------------------------------->

       CC  = ${CC}
    CFLAGS = ${CFLAGS}
       CXX = ${CXX}
  CXXFLAGS = ${CXXFLAGS}
  CPPFLAGS = ${CPPFLAGS}
   LDFLAGS = ${LDFLAGS}
  IFC_LIBS = ${IFC_LIBS}

    PREFIX = ${prefix}

<=============================================================>"
