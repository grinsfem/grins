BUILT_SOURCES   = .license.stamp
EXTRA_DIST =


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
toyprobdir = $(prefix)/examples/toyprob
toyprob_PROGRAMS = toy
toyprob_DATA = $(top_srcdir)/examples/toyprob/toy.in

toy_SOURCES = $(top_srcdir)/examples/toyprob/toy.C
toy_LDADD = $(LIBGRINS_LIBS)

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
#vortexexdir = $(prefix)/examples/thermal_vortex
#vortexex_PROGRAMS = vortex
#vortexex_DATA = $(top_srcdir)/examples/thermal_vortex/vortex.in

#vortex_SOURCES = $(top_srcdir)/examples/thermal_vortex/vortex.C
#vortex_LDADD = $(LIBGRINS_LIBS)

#=======================================================================
# Laminar Flame Example
#=======================================================================
laminarflamedir = $(prefix)/examples/laminar_flame
laminarflame_PROGRAMS = bunsen
laminarflame_DATA = $(top_srcdir)/examples/laminar_flame/bunsen.in \
                    $(top_srcdir)/examples/laminar_flame/elements.xml \
                    $(top_srcdir)/examples/laminar_flame/bunsen.xml

bunsen_SOURCES = $(top_srcdir)/examples/laminar_flame/bunsen.C \
                 $(top_srcdir)/examples/laminar_flame/constant_with_exp_layer.C \
                 $(top_srcdir)/examples/laminar_flame/bunsen_source.C \
                 $(top_srcdir)/examples/laminar_flame/ignite_initial_guess.C

bunsen_LDADD = $(LIBGRINS_LIBS)
bunsen_CPPFLAGS = -I$(top_srcdir)/examples/laminar_flame $(AM_CPPFLAGS)


all_sources = $(cavity_SOURCES)          \
              $(injection_SOURCES)       \
              $(rayleigh_taylor_SOURCES) \
              $(vortex_SOURCES)          \
              $(bunsen_SOURCES)



MAINTAINERCLEANFILES = Makefile.in

MOSTLYCLEANFILES = *.gcno

#---------------------------------
# Embedded license header support
#---------------------------------

.license.stamp: $(top_srcdir)/LICENSE
	$(top_srcdir)/src/common/lic_utils/update_license.pl $(top_srcdir)/LICENSE $(all_sources)
	echo 'updated source license headers' >$@
