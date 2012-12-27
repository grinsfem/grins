BUILT_SOURCES   = .license.stamp


#=======================================================================
# Cavity Benchmark
#=======================================================================
cavitybenchmarkdir = $(prefix)/examples/cavity_benchmark
cavitybenchmark_PROGRAMS = cavity
cavitybenchmark_DATA = $(top_srcdir)/examples/cavity_benchmark/cavity.in

cavity_SOURCES = $(top_srcdir)/examples/cavity_benchmark/cavity.C
cavity_LDADD = $(LIBGRINS_LIBS)


#=======================================================================
# Heated Channel Example
#=======================================================================


#=======================================================================
# Lid-Driven Cavity Example
#=======================================================================
liddrivencavitydir = $(prefix)/examples/lid_driven_cavity
liddrivencavity_DATA = $(top_srcdir)/examples/lid_driven_cavity/lid_driven_cavity.in

#=======================================================================
# Mass Injection Example
#=======================================================================
massinjectiondir = $(prefix)/examples/mass_injection
massinjection_PROGRAMS = injection
massinjection_DATA = $(top_srcdir)/examples/mass_injection/injection.in

injection_SOURCES = $(top_srcdir)/examples/mass_injection/injection.C
injection_LDADD = $(LIBGRINS_LIBS)

#=======================================================================
# Rayleigh-Taylor Example
#=======================================================================
rayleighdir = $(prefix)/examples/rayleigh_taylor
rayleigh_PROGRAMS = rayleigh
rayleigh_DATA = $(top_srcdir)/examples/rayleigh_taylor/rayleigh.in

rayleigh_SOURCES = $(top_srcdir)/examples/rayleigh_taylor/rayleigh.C
rayleigh_LDADD = $(LIBGRINS_LIBS)

#=======================================================================
# Reacting Low Mach Navier Stokes Example
#=======================================================================
rlmnsexdir = $(prefix)/examples/reacting_flow
rlmnsex_PROGRAMS = rlmns
rlmnsex_DATA = $(top_srcdir)/examples/reacting_flow/rlmns.in \
               $(top_srcdir)/examples/reacting_flow/elements.xml \
               $(top_srcdir)/examples/reacting_flow/air_2sp.xml

rlmns_SOURCES = $(top_srcdir)/examples/reacting_flow/rlmns.C
rlmns_LDADD = $(LIBGRINS_LIBS)

#=======================================================================
# Vortex Example
#=======================================================================
vortexexdir = $(prefix)/examples/thermal_vortex
vortexex_PROGRAMS = vortex
vortexex_DATA = $(top_srcdir)/examples/thermal_vortex/vortex.in

vortex_SOURCES = $(top_srcdir)/examples/thermal_vortex/vortex.C
vortex_LDADD = $(LIBGRINS_LIBS)

all_sources = $(cavity_SOURCES)          \
              $(injection_SOURCES)       \
              $(rayleigh_taylor_SOURCES) \
              $(vortex_SOURCES)

AM_CPPFLAGS = -I$(top_srcdir)/src/solver/include \
              -I$(top_srcdir)/src/physics/include \
              -I$(top_srcdir)/src/boundary_conditions/include \
              -I$(top_srcdir)/src/utilities/include \
              -I$(top_srcdir)/src/visualization/include \
              -I$(top_srcdir)/src/mat_props/include \
              -I$(top_srcdir)/src/bc_handling/include \
              -I$(top_srcdir)/src/qoi/include \
              -I$(top_srcdir)/src/kinetics/include

AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
               $(LIBMESH_CPPFLAGS) \
               $(GRVY_CFLAGS)

#----------------
# Cantera support
#----------------
if CANTERA_ENABLED
   AM_CPPFLAGS += $(CANTERA_CPPFLAGS)
   AM_LDFLAGS   = $(CANTERA_LDFLAGS)
endif

MAINTAINERCLEANFILES = Makefile.in

#---------------------------------
# Embedded license header support
#---------------------------------

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(all_sources)
	echo 'updated source license headers' >$@
