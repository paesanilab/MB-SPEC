include $(top_srcdir)/build-aux/programs.am

EXTRA_DIST = README.md

bin_PROGRAMS = mbspec-power-response mbspec-power-spectrum

mbspec_power_response_SOURCES = response.cpp
mbspec_power_response_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
mbspec_power_spectrum_SOURCES = spectrum.cpp
