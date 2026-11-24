include $(top_srcdir)/build-aux/common.am

LDADD = $(top_builddir)/src/tools/libtools.a $(MKL_LIBS)

noinst_PROGRAMS = response \
	       spectrum \
	       fresnel

EXTRA_DIST = \
    README.md \
    fresnel.json \
    fresnel_factors/refractive_index_water_k.txt \
    fresnel_factors/refractive_index_water_n.txt

response_SOURCES = response.cpp

spectrum_SOURCES = spectrum.cpp

fresnel_SOURCES = fresnel.cpp
fresnel_LDADD =
