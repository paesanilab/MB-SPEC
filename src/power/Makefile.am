include $(top_srcdir)/build-aux/common.am

LDADD = $(top_builddir)/src/tools/libtools.a $(LAPACKE_LIBS)

EXTRA_DIST = README.md

noinst_PROGRAMS = response \
               spectrum

response_SOURCES = response.cpp

spectrum_SOURCES = spectrum.cpp
