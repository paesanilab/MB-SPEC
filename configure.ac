AC_PREREQ([2.70])
AC_INIT([mb-spec], [1.0.0])
AC_CONFIG_SRCDIR([src/tools/io-tools.cpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign 1.11 no-define dist-bzip2 no-dist-gzip])
AM_SILENT_RULES([yes])

# Detect user-selected modes before AC_PROG_CXX can add its own standard flag.
AS_CASE([" $CXX $CPPFLAGS $CXXFLAGS "],
    [*" -std="* | *" /std:"*], [mb_default_cxx_standard=no],
    [mb_default_cxx_standard=yes])
AC_PROG_CXX
AC_LANG([C++])
AS_IF([test "x$mb_default_cxx_standard" = xyes],
    [CXXFLAGS="$CXXFLAGS -std=c++17"
     mb_required_cxx_standard=C++17
     mb_required_cxx_version=201703L],
    [mb_required_cxx_standard=C++11
     mb_required_cxx_version=201103L])
AC_MSG_CHECKING([for $mb_required_cxx_standard support])
AC_COMPILE_IFELSE(
    [AC_LANG_PROGRAM([[#include <array>]], [[
        static_assert(__cplusplus >= $mb_required_cxx_version, "the selected C++ standard is required");
        std::array<int, 1> values = {{1}};
        auto value = [&values]() { return values[0]; };
        return value();
    ]])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([MB-SPEC requires $mb_required_cxx_standard or newer in this configuration. Set CXX and CXXFLAGS accordingly.])])
AC_PROG_RANLIB

# OpenMP is required: the programs call its runtime functions directly.
AC_OPENMP
AS_IF([test "x$enable_openmp" = xno || test "x$ac_cv_prog_cxx_openmp" = xunsupported],
    [AC_MSG_ERROR([MB-SPEC requires a C++ compiler with OpenMP support; --disable-openmp is unsupported.])])

# Retain the legacy switch so disabling a required dependency fails clearly.
AC_ARG_ENABLE([mkl],
    [AS_HELP_STRING([--enable-mkl], [use required MKL LAPACKE routines (default=yes; disabling is unsupported)])],
    [], [enable_mkl=yes])
AS_IF([test "x$enable_mkl" != xyes],
    [AC_MSG_ERROR([MKL is required by MB-SPEC; --disable-mkl is unsupported.
Install Intel MKL and set --with-mkl-cflags and --with-mkl-libs for your installation.])])

AC_ARG_VAR([MKL_CFLAGS], [MKL compile flags])
AC_ARG_VAR([MKL_LIBS], [MKL linker flags and libraries])
AC_ARG_WITH([mkl-cflags],
    [AS_HELP_STRING([--with-mkl-cflags=FLAGS], [override MKL_CFLAGS])],
    [MKL_CFLAGS=$withval],
    [: ${MKL_CFLAGS='-m64 -I/opt/intel/mkl/include'}])
AC_ARG_WITH([mkl-libs],
    [AS_HELP_STRING([--with-mkl-libs=FLAGS], [override MKL_LIBS])],
    [MKL_LIBS=$withval],
    [: ${MKL_LIBS='-L/opt/intel/lib/intel64 -L/opt/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl'}])

saved_LIBS=$LIBS
saved_CXXFLAGS=$CXXFLAGS
LIBS="$MKL_LIBS $LIBS"
CXXFLAGS="$OPENMP_CXXFLAGS $MKL_CFLAGS $CXXFLAGS"
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
AS_IF([test "x$mkl_usable" != "xyes"], [
    AC_MSG_ERROR([MKL is required, but mkl_lapacke.h or the LAPACKE_dgetrf,
LAPACKE_dgetri, and LAPACKE_dsyev routines could not be compiled and linked.
Set --with-mkl-cflags and --with-mkl-libs for your Intel MKL installation.
See config.log for the compiler or linker error.])
])

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

AC_MSG_NOTICE([Configuration summary:
       CXX = $CXX
  CXXFLAGS = $CXXFLAGS
  CPPFLAGS = $CPPFLAGS
   LDFLAGS = $LDFLAGS
      LIBS = $LIBS
    OpenMP = $OPENMP_CXXFLAGS
MKL_CFLAGS = $MKL_CFLAGS
  MKL_LIBS = $MKL_LIBS
    prefix = $prefix])
