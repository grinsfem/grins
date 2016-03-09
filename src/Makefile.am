CLEANFILES =

#----------------------------------------
# Programs and libraries we want to build
#----------------------------------------

lib_LTLIBRARIES = libgrins.la

bin_PROGRAMS    = grins grins_version 

if CANTERA_ENABLED
   bin_PROGRAMS += cantera_kinetic_rates
endif

if ANTIOCH_ENABLED
   bin_PROGRAMS += antioch_thermo_tables
   bin_PROGRAMS += antioch_kinetic_rates
   bin_PROGRAMS += antioch_transport_values
endif

#----------------------------------------------
# List of source files to build dynamic library
#----------------------------------------------
libgrins_la_SOURCES =

# src/bc_handling files
libgrins_la_SOURCES += bc_handling/src/bc_handling_base.C
libgrins_la_SOURCES += bc_handling/src/generic_bc_handling.C
libgrins_la_SOURCES += bc_handling/src/heat_transfer_bc_handling.C
libgrins_la_SOURCES += bc_handling/src/axisym_heat_transfer_bc_handling.C
libgrins_la_SOURCES += bc_handling/src/inc_navier_stokes_bc_handling.C
libgrins_la_SOURCES += bc_handling/src/low_mach_navier_stokes_bc_handling.C
libgrins_la_SOURCES += bc_handling/src/reacting_low_mach_navier_stokes_bc_handling_instantiate.C
libgrins_la_SOURCES += bc_handling/src/solid_mechanics_bc_handling.C
libgrins_la_SOURCES += bc_handling/src/spalart_allmaras_bc_handling.C

# src/boundary_conditions files
libgrins_la_SOURCES += boundary_conditions/src/boundary_conditions.C
libgrins_la_SOURCES += boundary_conditions/src/bc_factory.C
libgrins_la_SOURCES += boundary_conditions/src/neumann_func_obj.C
libgrins_la_SOURCES += boundary_conditions/src/parabolic_profile.C
libgrins_la_SOURCES += boundary_conditions/src/gaussian_profile.C
libgrins_la_SOURCES += boundary_conditions/src/gaussian_xy_profile.C
libgrins_la_SOURCES += boundary_conditions/src/concentric_cylinder_profile.C
libgrins_la_SOURCES += boundary_conditions/src/dbc_container.C
libgrins_la_SOURCES += boundary_conditions/src/pressure_pinning.C
libgrins_la_SOURCES += boundary_conditions/src/pbc_container.C
libgrins_la_SOURCES += boundary_conditions/src/nbc_container.C
libgrins_la_SOURCES += boundary_conditions/src/catalytic_wall_instantiate.C
libgrins_la_SOURCES += boundary_conditions/src/catalycity_base.C
libgrins_la_SOURCES += boundary_conditions/src/constant_catalycity.C
libgrins_la_SOURCES += boundary_conditions/src/arrhenius_catalycity.C
libgrins_la_SOURCES += boundary_conditions/src/power_law_catalycity.C

# src/error_estimation files
libgrins_la_SOURCES += error_estimation/src/error_estimation_factory.C

# src/ic_handling files
libgrins_la_SOURCES += ic_handling/src/ic_handling_base.C
libgrins_la_SOURCES += ic_handling/src/generic_ic_handler.C

# src/physics files
libgrins_la_SOURCES += physics/src/multiphysics_sys.C
libgrins_la_SOURCES += physics/src/assembly_context.C
libgrins_la_SOURCES += physics/src/physics.C
libgrins_la_SOURCES += physics/src/stokes.C
libgrins_la_SOURCES += physics/src/inc_navier_stokes_base.C
libgrins_la_SOURCES += physics/src/inc_navier_stokes.C
libgrins_la_SOURCES += physics/src/inc_navier_stokes_stab_base.C
libgrins_la_SOURCES += physics/src/inc_navier_stokes_adjoint_stab.C
libgrins_la_SOURCES += physics/src/inc_navier_stokes_spgsm_stab.C
libgrins_la_SOURCES += physics/src/heat_transfer_base.C
libgrins_la_SOURCES += physics/src/heat_transfer.C
libgrins_la_SOURCES += physics/src/heat_transfer_source.C
libgrins_la_SOURCES += physics/src/heat_transfer_stab_base.C
libgrins_la_SOURCES += physics/src/heat_transfer_adjoint_stab.C
libgrins_la_SOURCES += physics/src/heat_transfer_spgsm_stab.C
libgrins_la_SOURCES += physics/src/axisym_heat_transfer.C
libgrins_la_SOURCES += physics/src/boussinesq_buoyancy_base.C
libgrins_la_SOURCES += physics/src/boussinesq_buoyancy.C
libgrins_la_SOURCES += physics/src/boussinesq_buoyancy_adjoint_stab.C
libgrins_la_SOURCES += physics/src/boussinesq_buoyancy_spgsm_stab.C
libgrins_la_SOURCES += physics/src/axisym_boussinesq_buoyancy.C
libgrins_la_SOURCES += physics/src/physics_factory.C
libgrins_la_SOURCES += physics/src/physics_factory_helper.C
libgrins_la_SOURCES += physics/src/low_mach_navier_stokes_base.C
libgrins_la_SOURCES += physics/src/low_mach_navier_stokes.C
libgrins_la_SOURCES += physics/src/low_mach_navier_stokes_stab_base.C
libgrins_la_SOURCES += physics/src/low_mach_navier_stokes_spgsm_stab.C
libgrins_la_SOURCES += physics/src/low_mach_navier_stokes_vms_stab.C
libgrins_la_SOURCES += physics/src/low_mach_navier_stokes_braack_stab.C
libgrins_la_SOURCES += physics/src/stab_helper.C
libgrins_la_SOURCES += physics/src/inc_navier_stokes_stab_helper.C
libgrins_la_SOURCES += physics/src/spalart_allmaras_stab_helper.C
libgrins_la_SOURCES += physics/src/heat_transfer_stab_helper.C
libgrins_la_SOURCES += physics/src/low_mach_navier_stokes_stab_helper.C
libgrins_la_SOURCES += physics/src/heat_conduction.C
libgrins_la_SOURCES += physics/src/reacting_low_mach_navier_stokes_abstract.C
libgrins_la_SOURCES += physics/src/reacting_low_mach_navier_stokes_instantiate.C
libgrins_la_SOURCES += physics/src/averaged_fan.C
libgrins_la_SOURCES += physics/src/averaged_fan_adjoint_stab.C
libgrins_la_SOURCES += physics/src/averaged_fan_base.C
libgrins_la_SOURCES += physics/src/averaged_turbine.C
libgrins_la_SOURCES += physics/src/averaged_turbine_adjoint_stab.C
libgrins_la_SOURCES += physics/src/averaged_turbine_base.C
libgrins_la_SOURCES += physics/src/scalar_ode.C
libgrins_la_SOURCES += physics/src/velocity_drag.C
libgrins_la_SOURCES += physics/src/velocity_drag_adjoint_stab.C
libgrins_la_SOURCES += physics/src/velocity_drag_base.C
libgrins_la_SOURCES += physics/src/velocity_penalty.C
libgrins_la_SOURCES += physics/src/velocity_penalty_adjoint_stab.C
libgrins_la_SOURCES += physics/src/velocity_penalty_base.C
libgrins_la_SOURCES += physics/src/parsed_velocity_source.C
libgrins_la_SOURCES += physics/src/parsed_velocity_source_adjoint_stab.C
libgrins_la_SOURCES += physics/src/parsed_velocity_source_base.C
libgrins_la_SOURCES += physics/src/turbulence_models_base.C
libgrins_la_SOURCES += physics/src/spalart_allmaras_parameters.C
libgrins_la_SOURCES += physics/src/spalart_allmaras_helper.C
libgrins_la_SOURCES += physics/src/spalart_allmaras.C
libgrins_la_SOURCES += physics/src/spalart_allmaras_stab_base.C
libgrins_la_SOURCES += physics/src/spalart_allmaras_spgsm_stab.C
libgrins_la_SOURCES += physics/src/solid_mechanics_abstract.C
libgrins_la_SOURCES += physics/src/elastic_membrane_abstract.C
libgrins_la_SOURCES += physics/src/elastic_membrane_base.C
libgrins_la_SOURCES += physics/src/elastic_membrane.C
libgrins_la_SOURCES += physics/src/elastic_membrane_rayleigh_damping.C
libgrins_la_SOURCES += physics/src/elastic_membrane_constant_pressure.C
libgrins_la_SOURCES += physics/src/elastic_cable_abstract.C
libgrins_la_SOURCES += physics/src/elastic_cable_base.C
libgrins_la_SOURCES += physics/src/elastic_cable.C
libgrins_la_SOURCES += physics/src/elastic_cable_rayleigh_damping.C
libgrins_la_SOURCES += physics/src/elastic_cable_constant_gravity.C
libgrins_la_SOURCES += physics/src/nonlinear_elasticity_instantiate.C
libgrins_la_SOURCES += physics/src/source_term_base.C
libgrins_la_SOURCES += physics/src/constant_source_term.C
libgrins_la_SOURCES += physics/src/parsed_source_term.C
libgrins_la_SOURCES += physics/src/physics_naming.C
libgrins_la_SOURCES += physics/src/convection_diffusion.C

# src/properties files
libgrins_la_SOURCES += properties/src/parsed_property_base.C
libgrins_la_SOURCES += properties/src/constant_conductivity.C
libgrins_la_SOURCES += properties/src/parsed_conductivity.C
libgrins_la_SOURCES += properties/src/constant_prandtl_conductivity.C
libgrins_la_SOURCES += properties/src/constant_specific_heat.C
libgrins_la_SOURCES += properties/src/constant_viscosity.C
libgrins_la_SOURCES += properties/src/parsed_viscosity.C
libgrins_la_SOURCES += properties/src/spalart_allmaras_viscosity.C
libgrins_la_SOURCES += properties/src/constant_source_func.C
libgrins_la_SOURCES += properties/src/cantera_mixture.C
libgrins_la_SOURCES += properties/src/cantera_thermo.C
libgrins_la_SOURCES += properties/src/cantera_transport.C
libgrins_la_SOURCES += properties/src/cantera_kinetics.C
libgrins_la_SOURCES += properties/src/cantera_evaluator.C
libgrins_la_SOURCES += properties/src/antioch_chemistry.C
libgrins_la_SOURCES += properties/src/antioch_mixture.C
libgrins_la_SOURCES += properties/src/antioch_kinetics.C
libgrins_la_SOURCES += properties/src/antioch_evaluator_instantiate.C
libgrins_la_SOURCES += properties/src/antioch_mixture_averaged_transport_mixture_instantiate.C
libgrins_la_SOURCES += properties/src/antioch_mixture_averaged_transport_evaluator_instantiate.C
libgrins_la_SOURCES += properties/src/antioch_constant_transport_instantiate.C
libgrins_la_SOURCES += properties/src/hookes_law.C
libgrins_la_SOURCES += properties/src/hookes_law_1d.C
libgrins_la_SOURCES += properties/src/hyperelasticity.C
libgrins_la_SOURCES += properties/src/incompressible_plane_stress_hyperelasticity.C
libgrins_la_SOURCES += properties/src/hyperelasticity_instantiate.C
libgrins_la_SOURCES += properties/src/mooney_rivlin.C
libgrins_la_SOURCES += properties/src/materials_parsing.C
libgrins_la_SOURCES += properties/src/viscosity_base.C

# src/qoi files
libgrins_la_SOURCES += qoi/src/average_nusselt_number.C
libgrins_la_SOURCES += qoi/src/qoi_base.C
libgrins_la_SOURCES += qoi/src/qoi_factory.C
libgrins_la_SOURCES += qoi/src/vorticity.C
libgrins_la_SOURCES += qoi/src/composite_qoi.C
libgrins_la_SOURCES += qoi/src/parsed_boundary_qoi.C
libgrins_la_SOURCES += qoi/src/parsed_interior_qoi.C
libgrins_la_SOURCES += qoi/src/weighted_flux_qoi.C

# src/solver files
libgrins_la_SOURCES += solver/src/grins_solver.C
libgrins_la_SOURCES += solver/src/mesh_builder.C
libgrins_la_SOURCES += solver/src/simulation.C
libgrins_la_SOURCES += solver/src/solver_factory.C
libgrins_la_SOURCES += solver/src/grins_steady_solver.C
libgrins_la_SOURCES += solver/src/grins_unsteady_solver.C
libgrins_la_SOURCES += solver/src/simulation_builder.C
libgrins_la_SOURCES += solver/src/solver_context.C
libgrins_la_SOURCES += solver/src/mesh_adaptive_solver_base.C
libgrins_la_SOURCES += solver/src/parameter_manager.C
libgrins_la_SOURCES += solver/src/parameter_user.C
libgrins_la_SOURCES += solver/src/steady_mesh_adaptive_solver.C
libgrins_la_SOURCES += solver/src/displacement_continuation_solver.C
libgrins_la_SOURCES += solver/src/solver_parsing.C
libgrins_la_SOURCES += solver/src/strategies_parsing.C
libgrins_la_SOURCES += solver/src/time_stepping_parsing.C
libgrins_la_SOURCES += solver/src/unsteady_mesh_adaptive_solver.C

#src/variables
libgrins_la_SOURCES += variables/src/variables_base.C
libgrins_la_SOURCES += variables/src/single_variable.C
libgrins_la_SOURCES += variables/src/fe_variables_base.C
libgrins_la_SOURCES += variables/src/single_fe_type_variable.C
libgrins_la_SOURCES += variables/src/velocity_variables.C
libgrins_la_SOURCES += variables/src/velocity_fe_variables.C
libgrins_la_SOURCES += variables/src/displacement_variables.C
libgrins_la_SOURCES += variables/src/displacement_fe_variables.C

# src/utilities files
libgrins_la_SOURCES += utilities/src/grins_version.C
libgrins_la_SOURCES += utilities/src/input_utils.C
libgrins_la_SOURCES += utilities/src/cached_values.C
libgrins_la_SOURCES += utilities/src/distance_function.C
libgrins_la_SOURCES += utilities/src/string_utils.C

# src/visualization files
libgrins_la_SOURCES += visualization/src/steady_visualization.C
libgrins_la_SOURCES += visualization/src/unsteady_visualization.C
libgrins_la_SOURCES += visualization/src/visualization.C
libgrins_la_SOURCES += visualization/src/visualization_factory.C
libgrins_la_SOURCES += visualization/src/postprocessed_quantities.C
libgrins_la_SOURCES += visualization/src/postprocessing_factory.C



# This is for installing .C files that are needed for instantiating
# with different template parameters in user code
sfilesdir = $(prefix)/src
sfiles_DATA  = physics/src/axisym_heat_transfer.C
sfiles_DATA += bc_handling/src/reacting_low_mach_navier_stokes_bc_handling.C
sfiles_DATA += boundary_conditions/src/catalytic_wall_base.C
sfiles_DATA += boundary_conditions/src/gas_recombination_catalytic_wall.C
sfiles_DATA += boundary_conditions/src/gas_solid_catalytic_wall.C
sfiles_DATA += physics/src/reacting_low_mach_navier_stokes.C
sfiles_DATA += properties/src/antioch_evaluator.C
sfiles_DATA += properties/src/antioch_mixture_averaged_transport_mixture.C
sfiles_DATA += properties/src/antioch_mixture_averaged_transport_evaluator.C
sfiles_DATA += properties/src/antioch_constant_transport_mixture.C
sfiles_DATA += properties/src/antioch_constant_transport_evaluator.C

#----------------------------
#INCLUDES we want distributed
#----------------------------
includedir = $(prefix)/include/grins
include_HEADERS =

# src/bc_handling headers
include_HEADERS += bc_handling/include/grins/axisym_heat_transfer_bc_handling.h
include_HEADERS += bc_handling/include/grins/bc_handling_base.h
include_HEADERS += bc_handling/include/grins/generic_bc_handling.h
include_HEADERS += bc_handling/include/grins/heat_transfer_bc_handling.h
include_HEADERS += bc_handling/include/grins/inc_navier_stokes_bc_handling.h
include_HEADERS += bc_handling/include/grins/low_mach_navier_stokes_bc_handling.h
include_HEADERS += bc_handling/include/grins/reacting_low_mach_navier_stokes_bc_handling.h
include_HEADERS += bc_handling/include/grins/solid_mechanics_bc_handling.h
include_HEADERS += bc_handling/include/grins/spalart_allmaras_bc_handling.h

# src/boundary_conditions headers
include_HEADERS += boundary_conditions/include/grins/bc_types.h
include_HEADERS += boundary_conditions/include/grins/boundary_conditions.h
include_HEADERS += boundary_conditions/include/grins/bc_factory.h
include_HEADERS += boundary_conditions/include/grins/neumann_func_obj.h
include_HEADERS += boundary_conditions/include/grins/parabolic_profile.h
include_HEADERS += boundary_conditions/include/grins/gaussian_profile.h
include_HEADERS += boundary_conditions/include/grins/gaussian_xy_profile.h
include_HEADERS += boundary_conditions/include/grins/concentric_cylinder_profile.h
include_HEADERS += boundary_conditions/include/grins/dbc_container.h
include_HEADERS += boundary_conditions/include/grins/pressure_pinning.h
include_HEADERS += boundary_conditions/include/grins/pbc_container.h
include_HEADERS += boundary_conditions/include/grins/nbc_container.h
include_HEADERS += boundary_conditions/include/grins/catalytic_wall_base.h
include_HEADERS += boundary_conditions/include/grins/gas_recombination_catalytic_wall.h
include_HEADERS += boundary_conditions/include/grins/gas_solid_catalytic_wall.h
include_HEADERS += boundary_conditions/include/grins/catalycity_base.h
include_HEADERS += boundary_conditions/include/grins/constant_catalycity.h
include_HEADERS += boundary_conditions/include/grins/arrhenius_catalycity.h
include_HEADERS += boundary_conditions/include/grins/power_law_catalycity.h

# src/common headers
include_HEADERS += common/include/grins/common.h

# src/error_estimation headers
include_HEADERS += error_estimation/include/grins/error_estimation_factory.h

# src/ic_handling headers
include_HEADERS += ic_handling/include/grins/ic_handling_base.h
include_HEADERS += ic_handling/include/grins/generic_ic_handler.h

# src/physics headers
include_HEADERS += physics/include/grins/multiphysics_sys.h
include_HEADERS += physics/include/grins/assembly_context.h
include_HEADERS += physics/include/grins/physics.h
include_HEADERS += physics/include/grins/var_typedefs.h
include_HEADERS += physics/include/grins/stokes.h
include_HEADERS += physics/include/grins/inc_navier_stokes_base.h
include_HEADERS += physics/include/grins/inc_navier_stokes.h
include_HEADERS += physics/include/grins/inc_navier_stokes_stab_base.h
include_HEADERS += physics/include/grins/inc_navier_stokes_adjoint_stab.h
include_HEADERS += physics/include/grins/inc_navier_stokes_spgsm_stab.h
include_HEADERS += physics/include/grins/heat_transfer_macros.h
include_HEADERS += physics/include/grins/heat_transfer_base.h
include_HEADERS += physics/include/grins/heat_transfer.h
include_HEADERS += physics/include/grins/heat_transfer_source.h
include_HEADERS += physics/include/grins/heat_transfer_stab_base.h
include_HEADERS += physics/include/grins/heat_transfer_adjoint_stab.h
include_HEADERS += physics/include/grins/heat_transfer_spgsm_stab.h
include_HEADERS += physics/include/grins/axisym_heat_transfer.h
include_HEADERS += physics/include/grins/boussinesq_buoyancy_base.h
include_HEADERS += physics/include/grins/boussinesq_buoyancy.h
include_HEADERS += physics/include/grins/boussinesq_buoyancy_adjoint_stab.h
include_HEADERS += physics/include/grins/boussinesq_buoyancy_spgsm_stab.h
include_HEADERS += physics/include/grins/axisym_boussinesq_buoyancy.h
include_HEADERS += physics/include/grins/physics_factory.h
include_HEADERS += physics/include/grins/physics_factory_helper.h
include_HEADERS += physics/include/grins/physics_naming.h
include_HEADERS += physics/include/grins/low_mach_navier_stokes_base.h
include_HEADERS += physics/include/grins/low_mach_navier_stokes.h
include_HEADERS += physics/include/grins/low_mach_navier_stokes_stab_base.h
include_HEADERS += physics/include/grins/low_mach_navier_stokes_spgsm_stab.h
include_HEADERS += physics/include/grins/low_mach_navier_stokes_vms_stab.h
include_HEADERS += physics/include/grins/low_mach_navier_stokes_braack_stab.h
include_HEADERS += physics/include/grins/stab_helper.h
include_HEADERS += physics/include/grins/inc_nav_stokes_macro.h
include_HEADERS += physics/include/grins/turbulence_models_macro.h
include_HEADERS += physics/include/grins/inc_navier_stokes_stab_helper.h
include_HEADERS += physics/include/grins/spalart_allmaras_stab_helper.h
include_HEADERS += physics/include/grins/heat_transfer_stab_helper.h
include_HEADERS += physics/include/grins/low_mach_navier_stokes_stab_helper.h
include_HEADERS += physics/include/grins/reacting_low_mach_navier_stokes_abstract.h
include_HEADERS += physics/include/grins/reacting_low_mach_navier_stokes_base.h
include_HEADERS += physics/include/grins/reacting_low_mach_navier_stokes.h
include_HEADERS += physics/include/grins/reacting_low_mach_navier_stokes_macro.h
include_HEADERS += physics/include/grins/heat_conduction.h
include_HEADERS += physics/include/grins/averaged_fan.h
include_HEADERS += physics/include/grins/averaged_fan_adjoint_stab.h
include_HEADERS += physics/include/grins/averaged_fan_base.h
include_HEADERS += physics/include/grins/averaged_turbine.h
include_HEADERS += physics/include/grins/averaged_turbine_adjoint_stab.h
include_HEADERS += physics/include/grins/averaged_turbine_base.h
include_HEADERS += physics/include/grins/scalar_ode.h
include_HEADERS += physics/include/grins/velocity_drag.h
include_HEADERS += physics/include/grins/velocity_drag_adjoint_stab.h
include_HEADERS += physics/include/grins/velocity_drag_base.h
include_HEADERS += physics/include/grins/velocity_penalty.h
include_HEADERS += physics/include/grins/velocity_penalty_adjoint_stab.h
include_HEADERS += physics/include/grins/velocity_penalty_base.h
include_HEADERS += physics/include/grins/parsed_velocity_source.h
include_HEADERS += physics/include/grins/parsed_velocity_source_adjoint_stab.h
include_HEADERS += physics/include/grins/parsed_velocity_source_base.h
include_HEADERS += physics/include/grins/turbulence_models_base.h
include_HEADERS += physics/include/grins/spalart_allmaras.h
include_HEADERS += physics/include/grins/spalart_allmaras_parameters.h
include_HEADERS += physics/include/grins/spalart_allmaras_helper.h
include_HEADERS += physics/include/grins/spalart_allmaras_stab_base.h
include_HEADERS += physics/include/grins/spalart_allmaras_spgsm_stab.h
include_HEADERS += physics/include/grins/solid_mechanics_abstract.h
include_HEADERS += physics/include/grins/elastic_membrane_abstract.h
include_HEADERS += physics/include/grins/elastic_membrane_base.h
include_HEADERS += physics/include/grins/elastic_membrane.h
include_HEADERS += physics/include/grins/elastic_cable_abstract.h
include_HEADERS += physics/include/grins/elastic_cable_base.h
include_HEADERS += physics/include/grins/elastic_cable.h
include_HEADERS += physics/include/grins/elastic_cable_rayleigh_damping.h
include_HEADERS += physics/include/grins/elastic_membrane_rayleigh_damping.h
include_HEADERS += physics/include/grins/elastic_membrane_constant_pressure.h
include_HEADERS += physics/include/grins/elastic_cable_constant_gravity.h
include_HEADERS += physics/include/grins/source_term_base.h
include_HEADERS += physics/include/grins/constant_source_term.h
include_HEADERS += physics/include/grins/parsed_source_term.h
include_HEADERS += physics/include/grins/antioch_instantiation_macro.h
include_HEADERS += physics/include/grins/convection_diffusion.h

# src/properties headers
include_HEADERS += properties/include/grins/property_types.h
include_HEADERS += properties/include/grins/parsed_property_base.h
include_HEADERS += properties/include/grins/constant_conductivity.h
include_HEADERS += properties/include/grins/parsed_conductivity.h
include_HEADERS += properties/include/grins/constant_prandtl_conductivity.h
include_HEADERS += properties/include/grins/constant_specific_heat.h
include_HEADERS += properties/include/grins/constant_viscosity.h
include_HEADERS += properties/include/grins/parsed_viscosity.h
include_HEADERS += properties/include/grins/spalart_allmaras_viscosity.h
include_HEADERS += properties/include/grins/turbulent_viscosity_macro.h
include_HEADERS += properties/include/grins/constant_source_func.h
include_HEADERS += properties/include/grins/cantera_mixture.h
include_HEADERS += properties/include/grins/cantera_thermo.h
include_HEADERS += properties/include/grins/cantera_transport.h
include_HEADERS += properties/include/grins/cantera_kinetics.h
include_HEADERS += properties/include/grins/cantera_evaluator.h
include_HEADERS += properties/include/grins/antioch_chemistry.h
include_HEADERS += properties/include/grins/antioch_kinetics.h
include_HEADERS += properties/include/grins/antioch_mixture.h
include_HEADERS += properties/include/grins/antioch_evaluator.h
include_HEADERS += properties/include/grins/antioch_mixture_averaged_transport_mixture.h
include_HEADERS += properties/include/grins/antioch_mixture_averaged_transport_evaluator.h
include_HEADERS += properties/include/grins/antioch_constant_transport_mixture.h
include_HEADERS += properties/include/grins/antioch_constant_transport_evaluator.h
include_HEADERS += properties/include/grins/elasticity_tensor.h
include_HEADERS += properties/include/grins/stress_strain_law.h
include_HEADERS += properties/include/grins/hookes_law.h
include_HEADERS += properties/include/grins/hookes_law_1d.h
include_HEADERS += properties/include/grins/hyperelasticity.h
include_HEADERS += properties/include/grins/incompressible_plane_stress_hyperelasticity.h
include_HEADERS += properties/include/grins/hyperelastic_strain_energy.h
include_HEADERS += properties/include/grins/mooney_rivlin.h
include_HEADERS += properties/include/grins/materials_parsing.h
include_HEADERS += properties/include/grins/viscosity_base.h

# src/qoi headers
include_HEADERS += qoi/include/grins/average_nusselt_number.h
include_HEADERS += qoi/include/grins/qoi_base.h
include_HEADERS += qoi/include/grins/qoi_factory.h
include_HEADERS += qoi/include/grins/qoi_names.h
include_HEADERS += qoi/include/grins/vorticity.h
include_HEADERS += qoi/include/grins/composite_qoi.h
include_HEADERS += qoi/include/grins/parsed_boundary_qoi.h
include_HEADERS += qoi/include/grins/parsed_interior_qoi.h
include_HEADERS += qoi/include/grins/weighted_flux_qoi.h

# src/solver headers
include_HEADERS += solver/include/grins/grins_solver.h
include_HEADERS += solver/include/grins/mesh_builder.h
include_HEADERS += solver/include/grins/simulation.h
include_HEADERS += solver/include/grins/solver_factory.h
include_HEADERS += solver/include/grins/grins_steady_solver.h
include_HEADERS += solver/include/grins/grins_unsteady_solver.h
include_HEADERS += solver/include/grins/simulation_builder.h
include_HEADERS += solver/include/grins/solver_context.h
include_HEADERS += solver/include/grins/mesh_adaptive_solver_base.h
include_HEADERS += solver/include/grins/parameter_manager.h
include_HEADERS += solver/include/grins/parameter_user.h
include_HEADERS += solver/include/grins/steady_mesh_adaptive_solver.h
include_HEADERS += solver/include/grins/displacement_continuation_solver.h
include_HEADERS += solver/include/grins/solver_parsing.h
include_HEADERS += solver/include/grins/solver_names.h
include_HEADERS += solver/include/grins/strategies_parsing.h
include_HEADERS += solver/include/grins/time_stepping_parsing.h
include_HEADERS += solver/include/grins/simulation_parsing.h
include_HEADERS += solver/include/grins/unsteady_mesh_adaptive_solver.h

#src/variables headers
include_HEADERS += variables/include/grins/single_variable.h
include_HEADERS += variables/include/grins/fe_variables_base.h
include_HEADERS += variables/include/grins/single_fe_type_variable.h
include_HEADERS += variables/include/grins/pressure_variable.h
include_HEADERS += variables/include/grins/pressure_fe_variable.h
include_HEADERS += variables/include/grins/velocity_variables.h
include_HEADERS += variables/include/grins/velocity_fe_variables.h
include_HEADERS += variables/include/grins/primitive_temp_variables.h
include_HEADERS += variables/include/grins/primitive_temp_fe_variables.h
include_HEADERS += variables/include/grins/displacement_variables.h
include_HEADERS += variables/include/grins/displacement_fe_variables.h
include_HEADERS += variables/include/grins/species_mass_fracs_variables.h
include_HEADERS += variables/include/grins/species_mass_fracs_fe_variables.h
include_HEADERS += variables/include/grins/thermo_pressure_variable.h
include_HEADERS += variables/include/grins/thermo_pressure_fe_variable.h
include_HEADERS += variables/include/grins/turbulence_variables.h
include_HEADERS += variables/include/grins/turbulence_fe_variables.h
include_HEADERS += variables/include/grins/variable_name_defaults.h
include_HEADERS += variables/include/grins/variables_base.h
include_HEADERS += variables/include/grins/variables_parsing.h

# src/utilities headers
include_HEADERS += $(top_builddir)/src/utilities/include/grins/grins_version.h
include_HEADERS += utilities/include/grins/grins_enums.h
include_HEADERS += utilities/include/grins/shared_ptr.h
include_HEADERS += utilities/include/grins/physical_constants.h
include_HEADERS += utilities/include/grins/input_utils.h
include_HEADERS += utilities/include/grins/math_constants.h
include_HEADERS += utilities/include/grins/cached_values.h
include_HEADERS += utilities/include/grins/cached_quantities_enum.h
include_HEADERS += utilities/include/grins/string_utils.h
include_HEADERS += utilities/include/grins/distance_function.h

# src/visualization headers
include_HEADERS += visualization/include/grins/steady_visualization.h
include_HEADERS += visualization/include/grins/unsteady_visualization.h
include_HEADERS += visualization/include/grins/visualization.h
include_HEADERS += visualization/include/grins/visualization_factory.h
include_HEADERS += visualization/include/grins/postprocessed_quantities.h
include_HEADERS += visualization/include/grins/postprocessing_factory.h

if LIBMESH_LIBTOOL
   libgrins_la_LIBADD = $(LIBMESH_LIBDIR)/libmesh_$(LIBMESH_METHOD).la
else
   libgrins_la_LIBADD = $(LIBMESH_LDFLAGS) $(LIBMESH_LIBS)
endif

if GRVY_ENABLED
   libgrins_la_LIBADD += $(GRVY_PREFIX)/lib/libgrvy.la
endif

if ANTIOCH_ENABLED
   libgrins_la_LIBADD += $(ANTIOCH_LDFLAGS) $(ANTIOCH_LIBS)
endif

#------------------------------------
# Source for stand alone GRINS driver
#------------------------------------

grins_SOURCES = solver/src/grins.C
grins_LDADD = libgrins.la
if !LIBMESH_LIBTOOL
grins_LDADD += $(LIBMESH_LDFLAGS) $(LIBMESH_LIBS)
endif

grins_version_SOURCES = utilities/src/version.C
grins_version_LDADD = libgrins.la
if !LIBMESH_LIBTOOL
grins_version_LDADD += $(LIBMESH_LDFLAGS) $(LIBMESH_LIBS)
endif

if CANTERA_ENABLED
   cantera_kinetic_rates_SOURCES = apps/cantera_kinetic_rates.C
   cantera_kinetic_rates_LDADD = libgrins.la
endif

if ANTIOCH_ENABLED
   antioch_thermo_tables_SOURCES = apps/antioch_thermo_tables.C
   antioch_thermo_tables_LDADD = libgrins.la

   antioch_kinetic_rates_SOURCES = apps/antioch_kinetic_rates.C
   antioch_kinetic_rates_LDADD = libgrins.la

   antioch_transport_values_SOURCES = apps/antioch_transport_values.C
   antioch_transport_values_LDADD = libgrins.la
endif

#--------------------------------------
#Local Directories to include for build
#--------------------------------------
AM_CPPFLAGS = 
AM_CPPFLAGS += -I$(top_srcdir)/src/bc_handling/include
AM_CPPFLAGS += -I$(top_srcdir)/src/boundary_conditions/include
AM_CPPFLAGS += -I$(top_srcdir)/src/common/include
AM_CPPFLAGS += -I$(top_srcdir)/src/error_estimation/include
AM_CPPFLAGS += -I$(top_srcdir)/src/ic_handling/include
AM_CPPFLAGS += -I$(top_srcdir)/src/physics/include
AM_CPPFLAGS += -I$(top_srcdir)/src/properties/include
AM_CPPFLAGS += -I$(top_srcdir)/src/qoi/include
AM_CPPFLAGS += -I$(top_srcdir)/src/solver/include
AM_CPPFLAGS += -I$(top_srcdir)/src/utilities/include
AM_CPPFLAGS += -I$(top_srcdir)/src/variables/include
AM_CPPFLAGS += -I$(top_srcdir)/src/visualization/include

AM_CPPFLAGS += -I$(top_builddir)/src/utilities/include #For grins_version.h

AM_CPPFLAGS += $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(LIBMESH_CPPFLAGS)
AM_CPPFLAGS += $(GRVY_CFLAGS)

AM_LDFLAGS =

#----------------
# Cantera support
#----------------
if CANTERA_ENABLED
   AM_CPPFLAGS += $(CANTERA_CPPFLAGS)
   AM_LDFLAGS += $(CANTERA_LDFLAGS)
   libgrins_la_LIBADD  += $(CANTERA_LIBS)
endif

#----------------
# Antioch support
#----------------
if ANTIOCH_ENABLED
   AM_CPPFLAGS += $(ANTIOCH_CPPFLAGS)
endif

#-------------
# MASA support
#-------------
if MASA_ENABLED
  AM_CPPFLAGS += $(MASA_CXXFLAGS)
  libgrins_la_LIBADD += $(MASA_LIBS)
endif

#--------------
# Code coverage
#--------------
if CODE_COVERAGE_ENABLED
  CLEANFILES += *.gcda *.gcno
endif

EXTRA_DIST  = common
EXTRA_DIST += $(sfiles_DATA)

# Required for AX_AM_MACROS
###@INC_AMINCLUDE@
