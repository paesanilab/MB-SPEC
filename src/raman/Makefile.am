SUBDIRS = .

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/tools $(MKL_CFLAGS)
AM_CXXFLAGS = $(CXX_FLAGS) -std=c++11

bin_PROGRAMS = polarizability \
	       response \
	       spectrum

EXTRA_DIST = README.md 

polarizability_SOURCES = polarizability.cpp 
polarizability_LDADD = ../tools/libtools.a $(MKL_LIBS)
polarizability_DEPS = ../tools/libtools.a

response_SOURCES = response.cpp 
response_LDADD = ../tools/libtools.a $(MKL_LIBS)
response_DEPS = ../tools/libtools.a

spectrum_SOURCES = spectrum.cpp 
spectrum_LDADD = ../tools/libtools.a $(MKL_LIBS)
spectrum_DEPS = ../tools/libtools.a
