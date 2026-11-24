AC_PREREQ([2.70])
AC_INIT([mb-spec], [1.0.0])
AC_CONFIG_SRCDIR([src/tools/io-tools.cpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign 1.16])
AM_SILENT_RULES([yes])

AC_PROG_CXX
# AM_PROG_AR also requests a C compiler for its archive-interface probe.
AM_PROG_AR
AC_PROG_RANLIB
AX_CXX_COMPILE_STDCXX([17], [noext], [mandatory])
AC_LANG([C++])

# OpenMP is required: the programs call its runtime functions directly.
AC_OPENMP
AS_IF([test "x$enable_openmp" = xno || test "x$ac_cv_prog_cxx_openmp" = xunsupported],
    [AC_MSG_ERROR([MB-SPEC requires a C++ compiler with OpenMP support; --disable-openmp is unsupported.])])

# Select the MKL installation and threading layer explicitly; LP64 is required.
AC_ARG_VAR([MKL_CFLAGS], [MKL compile flags])
AC_ARG_VAR([MKL_LIBS], [MKL linker flags and libraries])

saved_LIBS=$LIBS
saved_CXXFLAGS=$CXXFLAGS
LIBS="$MKL_LIBS $LIBS"
CXXFLAGS="$OPENMP_CXXFLAGS $MKL_CFLAGS $CXXFLAGS"
AC_MSG_CHECKING([for required MKL LAPACKE routines])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[#include <mkl_lapacke.h>]], [[
        double a[1] = {1.0};
        double eigenvalues[1];
        static_assert(sizeof(lapack_int) == sizeof(int), "LP64 MKL is required");
        int pivots[1];
        lapack_int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 1, 1, a, 1, pivots);
        info += LAPACKE_dgetri(LAPACK_ROW_MAJOR, 1, a, 1, pivots);
        info += LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', 1, a, 1, eigenvalues);
        return info;
    ]])],
    [AC_MSG_RESULT([yes])],
    [AC_MSG_RESULT([no])
     AC_MSG_ERROR([LP64 MKL is required, but mkl_lapacke.h or the LAPACKE_dgetrf,
LAPACKE_dgetri, and LAPACKE_dsyev routines could not be compiled and linked.
Set MKL_CFLAGS and MKL_LIBS for your Intel MKL installation and threading layer.
See config.log for the compiler or linker error.])])
CXXFLAGS=$saved_CXXFLAGS
LIBS=$saved_LIBS

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
