include $(top_srcdir)/build-aux/common.am

LDADD = $(top_builddir)/src/tools/libtools.a $(MKL_LIBS)

noinst_PROGRAMS = polarizability \
	       response \
	       spectrum

EXTRA_DIST = README.md

polarizability_SOURCES = polarizability.cpp

response_SOURCES = response.cpp

spectrum_SOURCES = spectrum.cpp
