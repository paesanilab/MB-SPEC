AM_CPPFLAGS = -I$(top_srcdir)/src/tools $(MKL_CFLAGS) \
              -DUNITTEST_DIR='"$(abs_srcdir)"' \
              -DUNITTEST_BUILDDIR='"$(abs_builddir)"' \
              -DUNITTEST_BINDIR='"$(abs_top_builddir)/src"'

AM_CXXFLAGS = -std=c++17 -fopenmp
LDADD = $(top_builddir)/src/tools/libtools.a $(MKL_LIBS) $(GSL_LIBS)

check_PROGRAMS = unittest-ir unittest-raman unittest-power unittest-sfg
TESTS = $(check_PROGRAMS)
AM_TESTS_ENVIRONMENT = MB_SPEC_EXTENDED=1; export MB_SPEC_EXTENDED;

unittest_ir_SOURCES = unittest-ir.cpp
unittest_raman_SOURCES = unittest-raman.cpp
unittest_power_SOURCES = unittest-power.cpp
unittest_sfg_SOURCES = unittest-sfg.cpp
noinst_HEADERS = test_support.h json_options.h

EXTRA_DIST = README.md generate_references.py reference_cases.py \
             inputs/ir.json inputs/raman.json inputs/power.json inputs/sfg.json inputs/fresnel.json \
             inputs/bulk-water.xyz inputs/bulk-water.lammpstrj \
             inputs/f-water.xyz inputs/f-water.lammpstrj \
             inputs/air-water.lammpstrj inputs/graphene-water.lammpstrj \
             inputs/alumina.xyz inputs/alumina.lammpstrj \
             references/ir.json references/raman.json references/power.json references/sfg.json
