include $(top_srcdir)/build-aux/common.am

LDADD = $(top_builddir)/src/tools/libtools.a $(MKL_LIBS)

EXTRA_DIST = README.md

bin_PROGRAMS = dipole \
	       response \
               spectrum

dipole_SOURCES = dipole.cpp

response_SOURCES = response.cpp

spectrum_SOURCES = spectrum.cpp
