include $(top_srcdir)/build-aux/common.am

noinst_LIBRARIES = libtools.a

check_PROGRAMS = test-mat-tools test-poly2b5
TESTS = $(check_PROGRAMS)
test_mat_tools_SOURCES = test-mat-tools.cpp
test_mat_tools_LDADD = libtools.a $(LAPACKE_LIBS)
test_poly2b5_SOURCES = test-poly2b5.cpp test-poly2b5-data.h
test_poly2b5_LDADD = libtools.a

EXTRA_DIST = \
    README.md \
    json/json.h \
    macros.h \
    struct-lammpstrj-frame.h \
    struct-sfg-pols.h

libtools_a_SOURCES = constants.h \
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
    poly2b5.h poly2b5.cpp poly2b5-internal.h \
    poly2b5-ha1.cpp poly2b5-ha2.cpp poly2b5-hb1.cpp poly2b5-hb2.cpp \
    poly2b5-oa.cpp poly2b5-ob.cpp \
    dip_poly2b5.h dip_poly2b5.cpp \
    lp-site.h lp-site.cpp \
    poly-l-2b3.h poly-l-2b3.cpp \
    alpha-el-poly2b3.h alpha-el-poly2b3.cpp \
    alpha-ttm4.h alpha-ttm4.cpp \
    system.h system.cpp
