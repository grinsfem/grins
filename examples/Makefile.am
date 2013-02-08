BUILT_SOURCES   = .license.stamp
EXTRA_DIST =

#=======================================================================
# Backward facing step
#=======================================================================
backwardfacingstepdir = $(prefix)/examples/backward_facing_step
backwardfacingstep_SCRIPTS = $(top_builddir)/examples/backward_facing_step/run.sh
backwardfacingstep_DATA  = $(top_srcdir)/examples/backward_facing_step/backward_facing_step.in
backwardfacingstep_DATA += $(top_srcdir)/examples/backward_facing_step/mesh.e

EXTRA_DIST += $(backwardfacingstep_DATA)

#=======================================================================
# Cavity Benchmark
#=======================================================================
cavitybenchmarkdir = $(prefix)/examples/cavity_benchmark
cavitybenchmark_PROGRAMS = cavity
cavitybenchmark_DATA  = $(top_srcdir)/examples/cavity_benchmark/cavity.in
cavitybenchmark_SCRIPTS = $(top_builddir)/examples/cavity_benchmark/run.sh

EXTRA_DIST += $(cavitybenchmark_DATA)


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
liddrivencavity_SCRIPTS = $(top_builddir)/examples/lid_driven_cavity/run.sh

EXTRA_DIST += $(liddrivencavity_DATA)

#=======================================================================
# Mass Injection Example
#=======================================================================
massinjectiondir = $(prefix)/examples/mass_injection
massinjection_PROGRAMS = injection
massinjection_DATA  = $(top_srcdir)/examples/mass_injection/injection.in
massinjection_DATA += $(top_srcdir)/examples/mass_injection/injection_grid.e
massinjection_SCRIPTS = $(top_builddir)/examples/mass_injection/run.sh

EXTRA_DIST += $(massinjection_DATA)

injection_SOURCES = $(top_srcdir)/examples/mass_injection/injection.C
injection_LDADD = $(LIBGRINS_LIBS)

#=======================================================================
# Rayleigh-Taylor Example
#=======================================================================
rayleighdir = $(prefix)/examples/rayleigh_taylor
rayleigh_PROGRAMS = rayleigh
rayleigh_DATA = $(top_srcdir)/examples/rayleigh_taylor/rayleigh.in
rayleigh_SCRIPTS = $(top_builddir)/examples/rayleigh_taylor/run.sh

EXTRA_DIST += $(rayleigh_DATA)

rayleigh_SOURCES = $(top_srcdir)/examples/rayleigh_taylor/rayleigh.C
rayleigh_LDADD = $(LIBGRINS_LIBS)

#=======================================================================
# Coupled Stokes-Navier-Stokes Example
#=======================================================================
stokesnsdir = $(prefix)/examples/coupled_stokes_navierstokes
stokesns_DATA  = $(top_srcdir)/examples/coupled_stokes_navierstokes/stokes_ns.in
stokesns_DATA += $(top_srcdir)/examples/coupled_stokes_navierstokes/ns.in
stokesns_DATA += $(top_srcdir)/examples/coupled_stokes_navierstokes/mesh.e
stokesns_SCRIPTS = $(top_builddir)/examples/coupled_stokes_navierstokes/run.sh

EXTRA_DIST += $(stokesns_DATA)

#=======================================================================
# Vortex Example
#=======================================================================
#vortexexdir = $(prefix)/examples/thermal_vortex
#vortexex_PROGRAMS = vortex
#vortexex_DATA = $(top_srcdir)/examples/thermal_vortex/vortex.in

#vortex_SOURCES = $(top_srcdir)/examples/thermal_vortex/vortex.C
#vortex_LDADD = $(LIBGRINS_LIBS)

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
              -I$(top_srcdir)/src/qoi/include

AM_CPPFLAGS += $(BOOST_CPPFLAGS) \
               $(LIBMESH_CPPFLAGS) \
               $(GRVY_CFLAGS)

MAINTAINERCLEANFILES = Makefile.in

MOSTLYCLEANFILES = *.gcno

#---------------------------------
# Embedded license header support
#---------------------------------

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(all_sources)
	echo 'updated source license headers' >$@
