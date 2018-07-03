//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-


// This class
#include "grins/multiphysics_sys.h"

// GRINS
#include "grins/assembly_context.h"
#include "grins/fe_variables_base.h"
#include "grins/variable_warehouse.h"
#include "grins/bc_builder.h"
#include "grins/constraint_builder.h"
#include "grins/composite_qoi.h"
#include "grins/nonlinear_solver_options.h"

// libMesh
#include "libmesh/composite_function.h"
#include "libmesh/getpot.h"
#include "libmesh/parameter_multiaccessor.h"

namespace GRINS
{

  MultiphysicsSystem::MultiphysicsSystem( libMesh::EquationSystems& es,
                                          const std::string& name,
                                          const unsigned int number )
    : FEMSystem(es, name, number),
      _use_numerical_jacobians_only(false)
  {}

  void MultiphysicsSystem::attach_physics_list( PhysicsList physics_list )
  {
    _physics_list = physics_list;
  }

  void MultiphysicsSystem::read_input_options( const GetPot& input )
  {
    // Cache this for building boundary condition later
    _input = &input;

    NonlinearSolverOptions nonlinear_solver_options(input);

    // Read options for MultiphysicsSystem first
    this->verify_analytic_jacobians = nonlinear_solver_options.verify_analytic_jacobians();
    this->print_solution_norms = input("screen-options/print_solution_norms", false );
    this->print_solutions = input("screen-options/print_solutions", false );
    this->print_residual_norms = input("screen-options/print_residual_norms", false );

    // backwards compatibility with old config files.
    /*! \todo Remove old print_residual nomenclature */
    this->print_residuals = input("screen-options/print_residual", false );
    if (this->print_residuals)
      libmesh_deprecated();

    this->print_residuals = input("screen-options/print_residuals", this->print_residuals );
    this->print_jacobian_norms = input("screen-options/print_jacobian_norms", false );
    this->print_jacobians = input("screen-options/print_jacobians", false );
    this->print_element_solutions = input("screen-options/print_element_solutions", false );
    this->print_element_residuals = input("screen-options/print_element_residuals", false );
    this->print_element_jacobians = input("screen-options/print_element_jacobians", false );

    _use_numerical_jacobians_only = nonlinear_solver_options.use_numerical_jacobians_only();

    // This is defined in libMesh::FEMSystem
    numerical_jacobian_h = nonlinear_solver_options.numerical_jacobian_h();

    nonlinear_solver_options.numerical_jacobian_h_vars_and_vals
      (_numerical_jacobian_h_variables, _numerical_jacobian_h_values);
  }

  void MultiphysicsSystem::init_data()
  {
    // Need this to be true because of our overloading of the
    // mass_residual function.
    // This is data in FEMSystem. MUST be set before FEMSystem::init_data.
    use_fixed_solution = true;

    // Initalize all the variables. We pass this pointer for the system.
    /* NOTE: We CANNOT fuse this loop with the others. This loop
       MUST complete first. */
    /*! \todo Figure out how to tell compilers not to fuse this loop when
      they want to be aggressive. */
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      {
        (physics_iter->second)->init_variables( this );
      }

    libmesh_assert(_input);
    BCBuilder::build_boundary_conditions(*_input,*this,_neumann_bcs);

    this->_constraint =
      ConstraintBuilder::build_constraint_object(*_input,*this);

    this->attach_constraint_object(*this->_constraint);

    // If any variables need custom numerical_jacobian_h, we can set those
    // values now that variable names are all registered with the System
    for (unsigned int i=0; i != _numerical_jacobian_h_values.size(); ++i)
      {
        unsigned int var_num =
          this->variable_number(_numerical_jacobian_h_variables[i]);
        this->set_numerical_jacobian_h_for_var
          (var_num, _numerical_jacobian_h_values[i]);
      }

    // Now set time_evolving variables
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      {
        (physics_iter->second)->set_time_evolving_vars( this );
      }

    // Set whether the problem we're solving is steady or not
    // Since the variable is static, just call one Physics class
    {
      PhysicsListIter physics_iter = _physics_list.begin();
      if (physics_iter != _physics_list.end() )
        (physics_iter->second)->set_is_steady((this->time_solver)->is_steady());
    }

    // Next, call parent init_data function to intialize everything.
    libMesh::FEMSystem::init_data();

    // After solution has been initialized we can project initial
    // conditions to it
    libMesh::CompositeFunction<libMesh::Number> ic_function;
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      {
        // Initialize builtin IC's for each physics
        (physics_iter->second)->init_ics( this, ic_function );
      }

    if (ic_function.n_subfunctions())
      {
        this->project_solution(&ic_function);
      }

    // Now do any auxillary initialization required by each Physics
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      {
        (physics_iter->second)->auxiliary_init( *this );
      }

    return;
  }

  std::unique_ptr<libMesh::DiffContext> MultiphysicsSystem::build_context()
  {
    AssemblyContext* context = new AssemblyContext(*this);

    std::unique_ptr<libMesh::DiffContext> ap(context);

    libMesh::DifferentiablePhysics* phys = libMesh::FEMSystem::get_physics();

    libmesh_assert(phys);

    // If we are solving a moving mesh problem, tell that to the Context
    context->set_mesh_system(phys->get_mesh_system());
    context->set_mesh_x_var(phys->get_mesh_x_var());
    context->set_mesh_y_var(phys->get_mesh_y_var());
    context->set_mesh_z_var(phys->get_mesh_z_var());

    ap->set_deltat_pointer( &deltat );

    // If we are solving the adjoint problem, tell that to the Context
    ap->is_adjoint() = this->get_time_solver().is_adjoint();

    return ap;
  }

  void MultiphysicsSystem::register_postprocessing_vars( const GetPot& input,
                                                         PostProcessedQuantities<libMesh::Real>& postprocessing )
  {
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      {
        (physics_iter->second)->register_postprocessing_vars( input, postprocessing );
      }

    return;
  }

  void MultiphysicsSystem::register_parameter
  ( const std::string & param_name,
    libMesh::ParameterMultiAccessor<libMesh::Number>& param_pointer )
  {
    //Loop over each physics to ask each for the requested parameter
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      {
        (physics_iter->second)->register_parameter( param_name,
                                                    param_pointer );
      }

    for (unsigned int bc=0; bc<_neumann_bcs.size(); bc++)
      _neumann_bcs[bc]->get_func()->register_parameter(param_name, param_pointer);
  }



  void MultiphysicsSystem::init_context( libMesh::DiffContext& context )
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    //Loop over each physics to initialize relevant variable structures for assembling system
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      {
        (physics_iter->second)->init_context( c );
      }

    return;
  }

  void MultiphysicsSystem::assembly( bool get_residual,
                                     bool get_jacobian,
                                     bool apply_heterogeneous_constraints,
                                     bool apply_no_constraints )
  {
    // First do any preassembly that the Physics requires (which by default is none)
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      (physics_iter->second)->preassembly(*this);

    // Now do the assembly
    libMesh::FEMSystem::assembly(get_residual,get_jacobian,
                                 apply_heterogeneous_constraints,
                                 apply_no_constraints);
  }

  void MultiphysicsSystem::reinit()
  {
    // First call Parent
    FEMSystem::reinit();

    // Now do per Physics reinit (which by default is none)
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      (physics_iter->second)->reinit(*this);

    // And now reinit the QoI
    if (this->qoi.size() > 0)
      {
        libMesh::DifferentiableQoI* diff_qoi = this->get_qoi();
        CompositeQoI* qoi = libMesh::cast_ptr<CompositeQoI*>(diff_qoi);
        qoi->reinit(*this);
      }
  }

  bool MultiphysicsSystem::_general_residual( bool request_jacobian,
                                              libMesh::DiffContext& context,
                                              ResFuncType resfunc,
                                              CacheFuncType cachefunc)
  {
    AssemblyContext& c = libMesh::cast_ref<AssemblyContext&>(context);

    bool compute_jacobian = true;
    if( !request_jacobian || _use_numerical_jacobians_only ) compute_jacobian = false;

    CachedValues & cache = c.get_cached_values();

    // Now compute cache for this element
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      {
        // shared_ptr gets confused by operator->*
        ((*(physics_iter->second)).*cachefunc)( c );
      }

    // Loop over each physics and compute their contributions
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      {
        if(c.has_elem())
          {
            if( (physics_iter->second)->enabled_on_elem( &c.get_elem() ) )
              {
                ((*(physics_iter->second)).*resfunc)( compute_jacobian, c );
              }
          }
        else
          {
            ((*(physics_iter->second)).*resfunc)( compute_jacobian, c );
          }
      }

    // We need to clear out the cache when we're done so we don't interfere
    // with other residual functions
    cache.clear();

    // TODO: Need to think about the implications of this because there might be some
    // TODO: jacobian terms we don't want to compute for efficiency reasons
    return compute_jacobian;
  }

  bool MultiphysicsSystem::element_time_derivative( bool request_jacobian,
                                                    libMesh::DiffContext& context )
  {
    return this->_general_residual
      (request_jacobian,
       context,
       &GRINS::Physics::element_time_derivative,
       &GRINS::Physics::compute_element_time_derivative_cache);
  }

  bool MultiphysicsSystem::side_time_derivative( bool request_jacobian,
                                                 libMesh::DiffContext& context )
  {
    bool jacobian_computed = this->apply_neumann_bcs(request_jacobian,
                                                     context);

    jacobian_computed = jacobian_computed &&
      this->_general_residual
      (request_jacobian,
       context,
       &GRINS::Physics::side_time_derivative,
       &GRINS::Physics::compute_side_time_derivative_cache);

    return jacobian_computed;
  }

  bool MultiphysicsSystem::nonlocal_time_derivative( bool request_jacobian,
                                                     libMesh::DiffContext& context )
  {
    return this->_general_residual
      (request_jacobian,
       context,
       &GRINS::Physics::nonlocal_time_derivative,
       &GRINS::Physics::compute_nonlocal_time_derivative_cache);
  }

  bool MultiphysicsSystem::element_constraint( bool request_jacobian,
                                               libMesh::DiffContext& context )
  {
    return this->_general_residual
      (request_jacobian,
       context,
       &GRINS::Physics::element_constraint,
       &GRINS::Physics::compute_element_constraint_cache);
  }

  bool MultiphysicsSystem::side_constraint( bool request_jacobian,
                                            libMesh::DiffContext& context )
  {
    return this->_general_residual
      (request_jacobian,
       context,
       &GRINS::Physics::side_constraint,
       &GRINS::Physics::compute_side_constraint_cache);
  }

  bool MultiphysicsSystem::nonlocal_constraint( bool request_jacobian,
                                                libMesh::DiffContext& context )
  {
    return this->_general_residual
      (request_jacobian,
       context,
       &GRINS::Physics::nonlocal_constraint,
       &GRINS::Physics::compute_nonlocal_constraint_cache);
  }

  bool MultiphysicsSystem::damping_residual( bool request_jacobian,
                                             libMesh::DiffContext& context )
  {
    return this->_general_residual
      (request_jacobian,
       context,
       &GRINS::Physics::damping_residual,
       &GRINS::Physics::compute_damping_residual_cache);
  }

  bool MultiphysicsSystem::mass_residual( bool request_jacobian,
                                          libMesh::DiffContext& context )
  {
    return this->_general_residual
      (request_jacobian,
       context,
       &GRINS::Physics::mass_residual,
       &GRINS::Physics::compute_mass_residual_cache);
  }

  bool MultiphysicsSystem::nonlocal_mass_residual( bool request_jacobian,
                                                   libMesh::DiffContext& context )
  {
    return this->_general_residual
      (request_jacobian,
       context,
       &GRINS::Physics::nonlocal_mass_residual,
       &GRINS::Physics::compute_nonlocal_mass_residual_cache);
  }

  std::shared_ptr<Physics> MultiphysicsSystem::get_physics( const std::string physics_name )
  {
    if( _physics_list.find( physics_name ) == _physics_list.end() )
      {
        std::cerr << "Error: Could not find physics " << physics_name << std::endl;
        libmesh_error();
      }

    return _physics_list[physics_name];
  }

  bool MultiphysicsSystem::has_physics( const std::string physics_name ) const
  {
    bool has_physics = false;

    if( _physics_list.find(physics_name) != _physics_list.end() )
      has_physics = true;

    return has_physics;
  }

  void MultiphysicsSystem::compute_postprocessed_quantity( unsigned int quantity_index,
                                                           const AssemblyContext& context,
                                                           const libMesh::Point& point,
                                                           libMesh::Real& value )
  {
    for( PhysicsListIter physics_iter = _physics_list.begin();
         physics_iter != _physics_list.end();
         physics_iter++ )
      {
        // Only compute if physics is active on current subdomain or globally
        if( (physics_iter->second)->enabled_on_elem( &context.get_elem() ) )
          {
            (physics_iter->second)->compute_postprocessed_quantity( quantity_index, context, point, value );
          }
      }
    return;
  }

  void MultiphysicsSystem::get_active_neumann_bcs( BoundaryID bc_id,
                                                   const std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs,
                                                   std::vector<std::shared_ptr<NeumannBCContainer> >& active_neumann_bcs )
  {
    // Manually writing the loop since std::copy_if is C++11 only
    for( std::vector<std::shared_ptr<NeumannBCContainer> >::const_iterator it = neumann_bcs.begin();
         it < neumann_bcs.end(); ++it )
      if( (*it)->has_bc_id( bc_id ) )
        active_neumann_bcs.push_back( *it );
  }

  bool MultiphysicsSystem::apply_neumann_bcs( bool request_jacobian,
                                              libMesh::DiffContext& context )
  {
    AssemblyContext& assembly_context =
      libMesh::cast_ref<AssemblyContext&>( context );

    std::vector<BoundaryID> ids;
    assembly_context.side_boundary_ids(ids);

    bool compute_jacobian = request_jacobian;
    if( !request_jacobian || _use_numerical_jacobians_only ) compute_jacobian = false;

    for( std::vector<BoundaryID>::const_iterator it = ids.begin();
         it != ids.end(); it++ )
      {
        BoundaryID bc_id = *it;

        libmesh_assert_not_equal_to(bc_id, libMesh::BoundaryInfo::invalid_id);

        // Retreive the NeumannBCContainers that are active on the current bc_id
        std::vector<std::shared_ptr<NeumannBCContainer> > active_neumann_bcs;
        this->get_active_neumann_bcs( bc_id, _neumann_bcs, active_neumann_bcs );

        if( !active_neumann_bcs.empty() )
          {
            typedef std::vector<std::shared_ptr<NeumannBCContainer> >::iterator BCIt;

            for( BCIt container = active_neumann_bcs.begin(); container < active_neumann_bcs.end(); ++container )
              {
                std::shared_ptr<NeumannBCAbstract>& func = (*container)->get_func();

                const FEVariablesBase& var = (*container)->get_fe_var();

                func->eval_flux( compute_jacobian, assembly_context,
                                 var.neumann_bc_sign(), Physics::is_axisymmetric() );
              }
          }
      } // end loop over boundary ids

    return compute_jacobian;
  }

} // namespace GRINS
