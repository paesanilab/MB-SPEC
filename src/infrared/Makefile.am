include $(top_srcdir)/build-aux/programs.am

EXTRA_DIST = README.md

bin_PROGRAMS = mbspec-ir-dipole mbspec-ir-response mbspec-ir-spectrum

mbspec_ir_dipole_SOURCES = dipole.cpp
mbspec_ir_dipole_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
# The LTP2011 1B dipole surface inverts 3x3 rotations through LAPACKE (mat_tools).
mbspec_ir_dipole_LDADD = $(LDADD) $(LAPACKE_LIBS)
mbspec_ir_response_SOURCES = response.cpp
mbspec_ir_spectrum_SOURCES = spectrum.cpp
