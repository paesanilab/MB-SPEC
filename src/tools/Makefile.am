include $(top_srcdir)/build-aux/common.am

# ps.cpp uses OpenMP threadprivate storage; mat_tools.cpp includes LAPACKE.
AM_CXXFLAGS = $(MBSPEC_CXXFLAGS) $(OPENMP_CXXFLAGS) $(LAPACKE_CFLAGS)

noinst_LIBRARIES = libtools.a libpoly.a

check_PROGRAMS = test-mat-tools test-poly2b5
TESTS = $(check_PROGRAMS)
test_mat_tools_SOURCES = test-mat-tools.cpp
test_mat_tools_LDADD = libtools.a $(LAPACKE_LIBS)
test_poly2b5_SOURCES = test-poly2b5.cpp test-poly2b5-data.h
test_poly2b5_LDADD = libpoly.a

EXTRA_DIST = \
    README.md \
    cli.h \
    json/json.h \
    macros.h \
    struct-lammpstrj-frame.h \
    struct-sfg-pols.h

libtools_a_SOURCES = constants.h \
    numeric_input.cpp numeric_input.h \
    spectral_smoothing.h \
    orient.cpp orient.h pbc-tools.cpp pbc-tools.h \
    mat_tools.cpp mat_tools.h \
    io-tools.cpp io-tools.h  \
    read-lammpstrj.cpp read-lammpstrj.h \
    ps.h ps.cpp gammq.h gammq.cpp \
    ltp2011.h ltp2011.cpp \
    ttm4-smear.h ttm4-smear.cpp ttm4-es.h ttm4-es.cpp \
    ttm4-periodic-es.h ttm4-periodic-es.cpp ttm4-dip.h ttm4-dip.cpp \
    avila.cpp avila.h \
    dip_poly2b5.h dip_poly2b5.cpp \
    lp-site.h lp-site.cpp \
    alpha-el-poly2b3.h alpha-el-poly2b3.cpp \
    alpha-ttm4.h alpha-ttm4.cpp \
    system.h system.cpp

# Machine-generated polynomial evaluators: thousands of lines of straight-line
# arithmetic on fixed-size arrays. They are built without --enable-sanitizers
# instrumentation, which makes them take tens of minutes to compile; their
# results are checked against the original polynomial by test-poly2b5 and by
# the workflow tests. Link libpoly.a after libtools.a.
libpoly_a_SOURCES = \
    poly2b5.h poly2b5.cpp poly2b5-internal.h \
    poly2b5-ha1.cpp poly2b5-ha2.cpp poly2b5-hb1.cpp poly2b5-hb2.cpp \
    poly2b5-oa.cpp poly2b5-ob.cpp \
    poly-l-2b3.h poly-l-2b3.cpp
libpoly_a_CXXFLAGS = $(WARN_CXXFLAGS)
