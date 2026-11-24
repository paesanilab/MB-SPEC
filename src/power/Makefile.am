include $(top_srcdir)/build-aux/common.am

LDADD = $(top_builddir)/src/tools/libtools.a $(MKL_LIBS)

EXTRA_DIST = README.md

bin_PROGRAMS = response \
               spectrum

response_SOURCES = response.cpp

spectrum_SOURCES = spectrum.cpp
