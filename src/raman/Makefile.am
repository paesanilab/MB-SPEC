include $(top_srcdir)/build-aux/programs.am

EXTRA_DIST = README.md

bin_PROGRAMS = mbspec-raman-polarizability mbspec-raman-response mbspec-raman-spectrum

mbspec_raman_polarizability_SOURCES = polarizability.cpp
mbspec_raman_polarizability_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
mbspec_raman_response_SOURCES = response.cpp
mbspec_raman_spectrum_SOURCES = spectrum.cpp
