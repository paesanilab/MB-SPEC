include $(top_srcdir)/build-aux/programs.am

AM_CPPFLAGS += -DUNITTEST_DIR='"$(abs_srcdir)"' \
               -DUNITTEST_BUILDDIR='"$(abs_builddir)"' \
               -DUNITTEST_BINDIR='"$(abs_top_builddir)/src"'

check_PROGRAMS = unittest-ir unittest-raman unittest-power unittest-sfg check-historical check-inputs
TESTS = $(check_PROGRAMS)
AM_TESTS_ENVIRONMENT = MB_SPEC_EXTENDED=1; export MB_SPEC_EXTENDED;

unittest_ir_SOURCES = unittest-ir.cpp
unittest_ir_LDADD = $(LDADD) $(LAPACKE_LIBS)
unittest_raman_SOURCES = unittest-raman.cpp
unittest_power_SOURCES = unittest-power.cpp
unittest_sfg_SOURCES = unittest-sfg.cpp
check_historical_SOURCES = check-historical.cpp
check_historical_LDADD =
check_inputs_SOURCES = check-inputs.cpp
check_inputs_LDADD =
noinst_HEADERS = test_support.h json_options.h

EXTRA_DIST = README.md generate_references.py reference_cases.py \
             inputs/ir.json inputs/raman.json inputs/power.json inputs/sfg.json inputs/fresnel.json \
             inputs/bulk-water.xyz inputs/bulk-water.lammpstrj \
             inputs/f-water.xyz inputs/f-water.lammpstrj \
             inputs/air-water.lammpstrj inputs/graphene-water.lammpstrj \
             inputs/alumina.xyz inputs/alumina.lammpstrj \
             references/ir.json references/raman.json references/power.json references/sfg.json \
             build-historical.sh historical

# Test runs retain diagnostics until the build tree is cleaned for distribution.
distclean-local:
	rm -rf .build/runs
