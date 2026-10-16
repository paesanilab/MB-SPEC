SUBDIRS = .

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/tools $(MKL_CFLAGS)

bin_PROGRAMS = response \
	       spectrum \
	       fresnel

EXTRA_DIST = \
    README.md \
    fresnel.json \
    fresnel_factors/refractive_index_water_k.txt \
    fresnel_factors/refractive_index_water_n.txt

response_SOURCES = response.cpp 
response_LDADD = ../tools/libtools.a $(MKL_LIBS)
response_DEPS = ../tools/libtools.a

spectrum_SOURCES = spectrum.cpp 
spectrum_LDADD = ../tools/libtools.a $(MKL_LIBS)
spectrum_DEPS = ../tools/libtools.a

fresnel_SOURCES = fresnel.cpp
