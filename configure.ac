AC_PREREQ([2.70])
AC_INIT([mb-spec], [1.0.0])
AC_PREFIX_DEFAULT([$(pwd -P)])
AC_CONFIG_SRCDIR([src/tools/io-tools.cpp])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([foreign 1.16])
AM_SILENT_RULES([yes])

AC_PROG_CXX
# AM_PROG_AR is deliberately omitted: it only supports non-POSIX archivers
# (Microsoft lib) and would require a C compiler for its probe.
AC_PROG_RANLIB
AX_CXX_COMPILE_STDCXX([17], [noext], [mandatory])
AC_LANG([C++])

# Enable common warnings for GCC-compatible compilers (GCC, Clang, icpx).
AC_ARG_ENABLE([werror],
    [AS_HELP_STRING([--enable-werror],
        [treat compiler warnings as errors @<:@default=no@:>@])],
    [], [enable_werror=no])
AS_IF([test "x$GXX" = xyes], [WARN_CXXFLAGS=-Wall], [WARN_CXXFLAGS=])
AS_IF([test "x$enable_werror" = xyes],
    [AS_IF([test "x$GXX" = xyes],
        [WARN_CXXFLAGS="$WARN_CXXFLAGS -Werror"],
        [AC_MSG_ERROR([--enable-werror requires a GCC-compatible C++ compiler.])])])
AC_SUBST([WARN_CXXFLAGS])

# Optional AddressSanitizer and UndefinedBehaviorSanitizer instrumentation.
AC_ARG_ENABLE([sanitizers],
    [AS_HELP_STRING([--enable-sanitizers],
        [build with AddressSanitizer and UndefinedBehaviorSanitizer @<:@default=no@:>@])],
    [], [enable_sanitizers=no])
SANITIZE_CXXFLAGS=
AS_IF([test "x$enable_sanitizers" = xyes],
    # UBSan only reports and continues by default; make it abort so that a
    # finding fails the tests.
    [AX_CHECK_LINK_FLAG([-fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer],
        [SANITIZE_CXXFLAGS="-fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer"],
        [AC_MSG_ERROR([--enable-sanitizers requires a compiler and linker supporting -fsanitize=address,undefined.])])])
AC_SUBST([SANITIZE_CXXFLAGS])

# OpenMP is required: the programs call its runtime functions directly.
AC_OPENMP
AS_IF([test "x$enable_openmp" = xno || test "x$ac_cv_prog_cxx_openmp" = xunsupported],
    [AC_MSG_ERROR([MB-SPEC requires a C++ compiler with OpenMP support; --disable-openmp is unsupported.])])

# LAPACKE is required; MKL is an explicitly selected implementation.
AC_ARG_WITH([lapack],
    [AS_HELP_STRING([--with-lapack=generic|mkl],
        [LAPACKE implementation to use @<:@default=generic@:>@])],
    [], [with_lapack=generic])
AS_CASE([$with_lapack],
    [generic], [],
    [mkl], [AC_DEFINE([USE_MKL], [1], [Use Intel MKL's LAPACKE header.])],
    [AC_MSG_ERROR([--with-lapack must be generic or mkl; LAPACKE cannot be disabled.])])

AC_ARG_VAR([LAPACKE_CFLAGS], [LAPACKE compile flags for the selected implementation])
AC_ARG_VAR([LAPACKE_LIBS], [LAPACKE linker flags and libraries for the selected implementation])
# Warn rather than fail: MKL environment modules may export these variables.
AS_IF([test "${MKL_CFLAGS+set}" = set || test "${MKL_LIBS+set}" = set],
    [AC_MSG_WARN([MKL_CFLAGS and MKL_LIBS are ignored. Use --with-lapack=mkl with LAPACKE_CFLAGS and LAPACKE_LIBS instead.])])

# pkg-config is optional; without it the pkg-config candidate is skipped.
# (pkg.m4's default action when pkg-config is missing is to abort.)
PKG_PROG_PKG_CONFIG([], [PKG_CONFIG=])
AS_IF([test "$with_lapack" = mkl],
    [lapacke_pc_module=mkl-dynamic-lp64-seq],
    [lapacke_pc_module=lapacke])

saved_LIBS=$LIBS
saved_CXXFLAGS=$CXXFLAGS
# An explicitly supplied LAPACKE_LIBS value, even empty, never falls back.
# Without one:
#   generic: pkg-config 'lapacke', then -llapacke, then -llapacke -llapack -lblas
#   mkl:     the compiler's own MKL support (e.g. -qmkl), then pkg-config
#            'mkl-dynamic-lp64-seq', then the sequential LP64 libraries under
#            $MKLROOT when it is set.
AS_IF([test "${LAPACKE_LIBS+set}" = set], [lapacke_candidates=user],
    [test "$with_lapack" = mkl && test -n "$MKLROOT"], [lapacke_candidates="user pkgconfig mklroot"],
    [test "$with_lapack" = mkl], [lapacke_candidates="user pkgconfig"],
    [lapacke_candidates="pkgconfig lapacke lapacke_with_blas"])
lapacke_found=no
for lapacke_candidate in $lapacke_candidates; do
    lapacke_try_CFLAGS=$LAPACKE_CFLAGS
    AS_IF([test "$lapacke_candidate" = pkgconfig],
        [AC_MSG_CHECKING([for pkg-config module $lapacke_pc_module])
         AS_IF([test -n "$PKG_CONFIG" && $PKG_CONFIG --exists "$lapacke_pc_module"],
             [AC_MSG_RESULT([yes])],
             [AC_MSG_RESULT([no])
              continue])])
    AS_CASE([$lapacke_candidate],
        [user], [lapacke_try_LIBS=$LAPACKE_LIBS],
        [pkgconfig],
        [AS_IF([test "${LAPACKE_CFLAGS+set}" != set],
             [lapacke_try_CFLAGS=`$PKG_CONFIG --cflags "$lapacke_pc_module"`])
         lapacke_try_LIBS=`$PKG_CONFIG --libs "$lapacke_pc_module"`],
        [mklroot],
        [AS_IF([test "${LAPACKE_CFLAGS+set}" != set],
             [lapacke_try_CFLAGS="-I${MKLROOT}/include"])
         lapacke_try_LIBS="-L${MKLROOT}/lib -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"],
        [lapacke], [lapacke_try_LIBS=-llapacke],
        [lapacke_with_blas], [lapacke_try_LIBS="-llapacke -llapack -lblas"])
    CXXFLAGS="$OPENMP_CXXFLAGS $lapacke_try_CFLAGS $saved_CXXFLAGS"
    LIBS="$lapacke_try_LIBS $saved_LIBS"
    AC_MSG_CHECKING([for LP64 $with_lapack LAPACKE with libraries '$lapacke_try_LIBS'])
    AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([[
#ifdef USE_MKL
#include <mkl_lapacke.h>
#else
#include <lapacke.h>
#endif
    ]], [[
        double a[1] = {1.0};
        double eigenvalues[1];
        static_assert(sizeof(lapack_int) == 4 && sizeof(lapack_int) == sizeof(int),
                      "LP64 LAPACKE with 32-bit integers is required");
        lapack_int pivots[1];
        lapack_int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, 1, 1, a, 1, pivots);
        info += LAPACKE_dgetri(LAPACK_ROW_MAJOR, 1, a, 1, pivots);
        info += LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', 1, a, 1, eigenvalues);
        return info;
    ]])],
    [AC_MSG_RESULT([yes])
     lapacke_found=yes
     LAPACKE_CFLAGS=$lapacke_try_CFLAGS
     LAPACKE_LIBS=$lapacke_try_LIBS
     break],
    [AC_MSG_RESULT([no])])
done
CXXFLAGS=$saved_CXXFLAGS
LIBS=$saved_LIBS
AS_IF([test "$lapacke_found" != yes],
    [AC_MSG_ERROR([LP64 $with_lapack LAPACKE is required, but its header or the
LAPACKE_dgetrf, LAPACKE_dgetri, and LAPACKE_dsyev routines could not be compiled
and linked. Set LAPACKE_CFLAGS and LAPACKE_LIBS for a matching LP64 installation.
Generic builds require lapacke.h; --with-lapack=mkl requires mkl_lapacke.h
(initialize your MKL environment so MKLROOT is set, or pass the flags).
Explicit LAPACKE_LIBS are used without fallback. See config.log for details.])])

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
  Warnings = $WARN_CXXFLAGS
Sanitizers = $SANITIZE_CXXFLAGS
    LAPACK = $with_lapack
LAPACKE_CFLAGS = $LAPACKE_CFLAGS
LAPACKE_LIBS = $LAPACKE_LIBS
    prefix = $prefix])
