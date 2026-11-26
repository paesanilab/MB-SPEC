include $(top_srcdir)/build-aux/programs.am

EXTRA_DIST = README.md

bin_PROGRAMS = mbspec-sfg-response mbspec-sfg-spectrum mbspec-sfg-fresnel

mbspec_sfg_response_SOURCES = response.cpp
mbspec_sfg_response_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
mbspec_sfg_spectrum_SOURCES = spectrum.cpp

# fresnel is self-contained and does not need libtools or LAPACKE.
mbspec_sfg_fresnel_SOURCES = fresnel.cpp
mbspec_sfg_fresnel_LDADD =

# Example Fresnel configuration and water refractive-index tables, installed
# to $(pkgdatadir)/sfg/.
sfgdatadir = $(pkgdatadir)/sfg
nobase_dist_sfgdata_DATA = \
    fresnel.json \
    fresnel_factors/refractive_index_water_k.txt \
    fresnel_factors/refractive_index_water_n.txt
