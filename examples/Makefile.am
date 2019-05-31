EXTRA_DIST =

AM_CPPFLAGS =
AM_CPPFLAGS += -I$(top_srcdir)/src/boundary_conditions/include
AM_CPPFLAGS += -I$(top_srcdir)/src/common/include
AM_CPPFLAGS += -I$(top_srcdir)/src/physics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/properties/include
AM_CPPFLAGS += -I$(top_srcdir)/src/qoi/include
AM_CPPFLAGS += -I$(top_srcdir)/src/solver/include
AM_CPPFLAGS += -I$(top_srcdir)/src/strategies/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_srcdir)/src/variables/include
AM_CPPFLAGS += -I$(top_srcdir)/src/visualization/include

AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)

AM_LDFLAGS =
#----------------
# Cantera support
#----------------
if CANTERA_ENABLED
   AM_CPPFLAGS += $(CANTERA_CPPFLAGS)
   AM_LDFLAGS  += $(CANTERA_LDFLAGS)
   LIBGRINS_LIBS += $(CANTERA_LIBS)
endif

#----------------
# Antioch support
#----------------
if ANTIOCH_ENABLED
   AM_CPPFLAGS += $(ANTIOCH_CPPFLAGS)
endif

#=======================================================================
# Backward facing step
#=======================================================================
backwardfacingstepdir = $(prefix)/examples/backward_facing_step
backwardfacingstep_SCRIPTS = $(top_builddir)/examples/backward_facing_step/run.sh
backwardfacingstep_DATA  = $(top_srcdir)/examples/backward_facing_step/backward_facing_step.in
backwardfacingstep_DATA += $(top_srcdir)/examples/backward_facing_step/mesh.e
backwardfacingstep_DATA += $(top_srcdir)/examples/backward_facing_step/README

EXTRA_DIST += $(backwardfacingstep_DATA)

#=======================================================================
# SA Airfoil
#=======================================================================
saairfoildir = $(prefix)/examples/sa_airfoil
saairfoil_SCRIPTS = $(top_builddir)/examples/sa_airfoil/run.sh
saairfoil_DATA  = $(top_srcdir)/examples/sa_airfoil/sa_airfoil.in
saairfoil_DATA += $(top_srcdir)/examples/sa_airfoil/naca_0012_medium_boundary_id.xda
saairfoil_DATA += $(top_srcdir)/examples/sa_airfoil/turbulent_airfoil_Re_6000000.xda

EXTRA_DIST += $(saairfoil_DATA)

#=======================================================================
# Cavity Benchmark
#=======================================================================
cavitybenchmarkdir = $(prefix)/examples/cavity_benchmark
cavitybenchmark_PROGRAMS = cavity
cavitybenchmark_DATA  = $(top_srcdir)/examples/cavity_benchmark/cavity.in
cavitybenchmark_SCRIPTS = $(top_builddir)/examples/cavity_benchmark/run.sh

EXTRA_DIST += $(cavitybenchmark_DATA)

cavity_SOURCES = cavity_benchmark/cavity.C
cavity_LDADD = $(LIBGRINS_LIBS)


#=======================================================================
# Lid-Driven Cavity Example
#=======================================================================
liddrivencavitydir = $(prefix)/examples/lid_driven_cavity
liddrivencavity_DATA = $(top_srcdir)/examples/lid_driven_cavity/lid_driven_cavity.in
liddrivencavity_SCRIPTS = $(top_builddir)/examples/lid_driven_cavity/run.sh

EXTRA_DIST += $(liddrivencavity_DATA)

#=======================================================================
# Velocity Penalty Example
#=======================================================================
velocitypenaltydir = $(prefix)/examples/velocity_penalty
velocitypenalty_DATA = $(top_srcdir)/examples/velocity_penalty/velocity_penalty.in
velocitypenalty_SCRIPTS = $(top_builddir)/examples/velocity_penalty/run.sh

EXTRA_DIST += $(velocitypenalty_DATA)

#=======================================================================
# Convection Cell Example
#=======================================================================
convectioncelldir = $(prefix)/examples/convection_cell
convectioncell_DATA = $(top_srcdir)/examples/convection_cell/convection_cell.in
convectioncell_SCRIPTS = $(top_builddir)/examples/convection_cell/run.sh

EXTRA_DIST += $(convectioncell_DATA)

#=======================================================================
# Rayleigh-Taylor Example
#=======================================================================
rayleighdir = $(prefix)/examples/rayleigh_taylor
rayleigh_DATA = $(top_srcdir)/examples/rayleigh_taylor/rayleigh.in
rayleigh_SCRIPTS = $(top_builddir)/examples/rayleigh_taylor/run.sh

EXTRA_DIST += $(rayleigh_DATA)

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
# Averaged Fan/Turbine Examples
#=======================================================================
simplefandir = $(prefix)/examples/simple_fan
simplefan_DATA = $(top_srcdir)/examples/simple_fan/fan.in
simplefan_SCRIPTS = $(top_builddir)/examples/simple_fan/run.sh

EXTRA_DIST += $(simplefan_DATA)

simplepropdir = $(prefix)/examples/simple_prop
simpleprop_DATA = $(top_srcdir)/examples/simple_prop/fan.in
simpleprop_SCRIPTS = $(top_builddir)/examples/simple_prop/run.sh

EXTRA_DIST += $(simpleprop_DATA)


#=======================================================================
# Elastic Sheet Example
#=======================================================================
elasticsheetdir = $(prefix)/examples/elastic_sheet
elasticsheet_PROGRAMS = stretching_sheet
elasticsheet_DATA  = $(top_srcdir)/examples/elastic_sheet/sheet.in
elasticsheet_DATA += $(top_srcdir)/examples/elastic_sheet/displacement_continuation_solver.h
elasticsheet_SCRIPTS = $(top_builddir)/examples/elastic_sheet/run.sh

EXTRA_DIST += $(elasticsheet_DATA)

stretching_sheet_SOURCES  = elastic_sheet/stretching_sheet.C
stretching_sheet_SOURCES += elastic_sheet/displacement_continuation_solver.C

stretching_sheet_LDADD = $(LIBGRINS_LIBS)
stretching_sheet_CPPFLAGS = -I$(top_srcdir)/examples/elastic_sheet $(AM_CPPFLAGS)

#=======================================================================
# Inflating Sheet Example
#=======================================================================
inflatingsheetdir = $(prefix)/examples/inflating_sheet
inflatingsheet_PROGRAMS = sheet
inflatingsheet_DATA  = $(top_srcdir)/examples/inflating_sheet/sheet.in
inflatingsheet_DATA += $(top_srcdir)/examples/inflating_sheet/sheet_coarse.e
inflatingsheet_DATA += $(top_srcdir)/examples/inflating_sheet/pressure_continuation_solver.h
inflatingsheet_SCRIPTS = $(top_builddir)/examples/inflating_sheet/run.sh

EXTRA_DIST += $(inflatingsheet_DATA)

sheet_SOURCES  = inflating_sheet/inflating_sheet.C
sheet_SOURCES += inflating_sheet/pressure_continuation_solver.C

sheet_LDADD = $(LIBGRINS_LIBS)
sheet_CPPFLAGS = -I$(top_srcdir)/examples/inflating_sheet $(AM_CPPFLAGS)

#=======================================================================
# Elastic Cable Example
#=======================================================================
cabledispdir = $(prefix)/examples/cable_displacement
cabledisp_DATA = $(top_srcdir)/examples/cable_displacement/cable.in
cabledisp_SCRIPTS = $(top_builddir)/examples/cable_displacement/run.sh

EXTRA_DIST += $(cabledisp_DATA)

#=======================================================================
# Suspended Cable Example
#=======================================================================
suspendedcabledir = $(prefix)/examples/suspended_cable
suspendedcable_DATA  = $(top_srcdir)/examples/suspended_cable/cable.in
suspendedcable_SCRIPTS = $(top_builddir)/examples/suspended_cable/run.sh

EXTRA_DIST += $(suspendedcable_DATA)

#=======================================================================
# Inflating Quarter Sheet with 1D Stiffeners Example
#=======================================================================
mixeddiminflatingsheetdir = $(prefix)/examples/mixed_dim_inflating_sheet
mixeddiminflatingsheet_DATA  = $(top_srcdir)/examples/mixed_dim_inflating_sheet/sheet.in
mixeddiminflatingsheet_DATA += $(top_srcdir)/examples/mixed_dim_inflating_sheet/quarter_sheet.exo
mixeddiminflatingsheet_DATA += $(top_srcdir)/examples/mixed_dim_inflating_sheet/README
mixeddiminflatingsheet_SCRIPTS = $(top_builddir)/examples/mixed_dim_inflating_sheet/run.sh

EXTRA_DIST += $(mixeddiminflatingsheet_DATA)


#=======================================================================
# Transient AMR Example
#=======================================================================
transientamrdir = $(prefix)/examples/transient_amr
transientamr_DATA  = $(top_srcdir)/examples/transient_amr/convection_diffusion_unsteady_2d_amr.in
transientamr_DATA += $(top_srcdir)/examples/transient_amr/mixed_quad_tri_square_mesh.xda
transientamr_SCRIPTS = $(top_builddir)/examples/transient_amr/run.sh

EXTRA_DIST += $(transientamr_DATA)

#=======================================================================
# Ozone Flame Example
#=======================================================================
ozonedir = $(prefix)/examples/ozone_flame
ozone_DATA  = $(top_srcdir)/examples/ozone_flame/common_input.in
ozone_DATA += $(top_srcdir)/examples/ozone_flame/ozone_antioch.in
ozone_DATA += $(top_srcdir)/examples/ozone_flame/ozone_cantera.in
ozone_DATA += $(top_srcdir)/examples/ozone_flame/cantera_steady.in
ozone_DATA += $(top_srcdir)/examples/ozone_flame/ozone.xml
ozone_DATA += $(top_srcdir)/examples/ozone_flame/ozone_cea_data.dat
ozone_DATA += $(top_srcdir)/examples/ozone_flame/ozone_species_data.dat
ozone_DATA += $(top_srcdir)/examples/ozone_flame/restart.xdr
ozone_DATA += $(top_srcdir)/examples/ozone_flame/restart_steady.in
ozone_DATA += $(top_srcdir)/examples/ozone_flame/restart_amr.in
ozone_DATA += $(top_srcdir)/examples/ozone_flame/README
ozone_SCRIPTS  = $(top_builddir)/examples/ozone_flame/run.sh
ozone_SCRIPTS += $(top_builddir)/examples/ozone_flame/crun.sh

EXTRA_DIST += $(ozone_DATA)

#=======================================================================
# Spectroscopy Example
#=======================================================================
spectroscopydir = $(prefix)/examples/spectroscopy
spectroscopy_PROGRAMS = spec
spectroscopy_DATA  = $(top_srcdir)/examples/spectroscopy/spectroscopy.in
spectroscopy_DATA += $(top_srcdir)/examples/spectroscopy/opt_spectroscopy.in
spectroscopy_DATA += $(top_srcdir)/examples/spectroscopy/CO_data.dat
spectroscopy_DATA += $(top_srcdir)/examples/spectroscopy/CO_O_O2.xml
spectroscopy_DATA += $(top_srcdir)/examples/spectroscopy/CO_partition_function.dat
spectroscopy_DATA += $(top_srcdir)/examples/spectroscopy/README
spectroscopy_SCRIPTS = $(top_builddir)/examples/spectroscopy/run.sh

EXTRA_DIST += $(spectroscopy_DATA)

spec_SOURCES  = spectroscopy/spectroscopy.C

spec_LDADD = $(LIBGRINS_LIBS)

#=======================================================================
# Multigrid Examples
#=======================================================================
multigriddir       = $(prefix)/examples/multigrid
multigrid_DATA     = $(top_srcdir)/examples/multigrid/forced_poisson.in
multigrid_DATA    += $(top_srcdir)/examples/multigrid/stokes_fieldsplit.in
multigrid_DATA    += $(top_srcdir)/examples/multigrid/thermo_coupled_flow.in
multigrid_DATA    += $(top_srcdir)/examples/multigrid/README
multigrid_SCRIPTS  = $(top_builddir)/examples/multigrid/run_poisson.sh
multigrid_SCRIPTS += $(top_builddir)/examples/multigrid/run_stokes.sh
multigrid_SCRIPTS += $(top_builddir)/examples/multigrid/run_thermo_coupled_flow.sh

EXTRA_DIST += $(multigrid_DATA)

#=======================================================================
# Torsion Hyperelasticity Example
#=======================================================================
torsionhyperelasticitydir = $(prefix)/examples/torsion_hyperelasticity
torsionhyperelasticity_SCRIPTS = $(top_builddir)/examples/torsion_hyperelasticity/run.sh
torsionhyperelasticity_DATA  = $(top_srcdir)/examples/torsion_hyperelasticity/torsion.in
torsionhyperelasticity_DATA += $(top_srcdir)/examples/torsion_hyperelasticity/inc_torsion.in

EXTRA_DIST += $(torsionhyperelasticity_DATA)

MAINTAINERCLEANFILES = Makefile.in

MOSTLYCLEANFILES = *.gcno
