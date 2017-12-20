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


#ifndef GRINS_MULTIPHYSICS_SYS_H
#define GRINS_MULTIPHYSICS_SYS_H

// C++
#include <string>

// GRINS
#include "grins_config.h"
#include "grins/physics.h"
#include "grins/neumann_bc_container.h"

// libMesh
#include "libmesh/fem_system.h"

// libMesh forward declartions
class GetPot;

namespace libMesh
{
  class EquationSystems;
  class DiffContext;

  template <typename Scalar>
  class ParameterMultiAccessor;
}

namespace GRINS
{
  // Forward Declarations
  template <typename Scalar>
  class PostProcessedQuantities;

  //! Interface with libMesh for solving Multiphysics problems.
  /*!
    MultiphysicsSystem (through libMesh::FEMSystem) solves the following equation:

    \f$M(u)\dot{u} = F(u)\f$

    M = mass matrix
    u = solution vector
    F = time derivative

    Note that for the nonlinear system that is solved for implicit
    time stepping is:

    \f$M(u_{\theta})(u^n - u^{n+1}) + \Delta t F(u) = 0\f$

    *_time_derivative correspond to calculating terms for \f$F(u)\f$
    *_mass_residual correspond to calculating terms for \f$M(u)\dot{u}\f$
    */
  //TODO: is it F(u) or F(u_{\theta})?
  class MultiphysicsSystem : public libMesh::FEMSystem
  {
  public:

    //! Constructor. Will be called by libMesh only.
    MultiphysicsSystem( libMesh::EquationSystems& es,
                        const std::string& name,
                        const unsigned int number );

    //! Destructor. Clean up all physics allocations.
    ~MultiphysicsSystem(){};

    //! PhysicsList gets built by GRINS::PhysicsFactory and attached here.
    void attach_physics_list( PhysicsList physics_list );

    const PhysicsList& get_physics_list() const
    { return _physics_list; }

    //! Reads input options for this class and all physics that are enabled
    /*!
      This function reads the input options for the MultiphysicsSystem class and then
      enables each of the requested physics in the system. Finally, the input options
      for each of the physics will be read.
    */
    virtual void read_input_options( const GetPot& input );

    //! System initialization. Calls each physics implementation of init_variables()
    virtual void init_data();

    //! Each Physics will register their postprocessed quantities with this call
    void register_postprocessing_vars( const GetPot& input,
                                       PostProcessedQuantities<libMesh::Real>& postprocessing );

    //! Each Physics will register its copy(s) of an independent variable
    //  named in this call.
    void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number>& param_pointer );

    //! Override FEMSystem::build_context in order to use our own AssemblyContext
    virtual std::unique_ptr<libMesh::DiffContext> build_context();

    //! Context initialization. Calls each physics implementation of init_context()
    virtual void init_context( libMesh::DiffContext &context );

    //! Override FEMSystem::assembly
    /*! This allows us to insert things like preassembly(). */
    virtual void assembly( bool get_residual,
                           bool get_jacobian,
                           bool apply_heterogeneous_constraints = false,
                           bool apply_no_constraints = false );

    //! Override FEMSystem::reinit
    /*! This will allow each Physics to reinit things internally that need it,
      such as point locators. */
    virtual void reinit();

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    //! Element interior contributions to \f$F(u)\f$ which have time varying components.
    virtual bool element_time_derivative( bool request_jacobian,
                                          libMesh::DiffContext& context );

    //! Boundary contributions to \f$F(u)\f$ which have time varying components.
    virtual bool side_time_derivative( bool request_jacobian,
                                       libMesh::DiffContext& context );

    //! Contributions to \f$F(u)\f$ on SCALAR variables which have time varying components.
    virtual bool nonlocal_time_derivative( bool request_jacobian,
                                           libMesh::DiffContext& context );

    //! Element interior contributions to \f$F(u)\f$ which do not have time varying components.
    //! Element interior contributions to \f$F(u)\f$ which do not have time varying components.
    virtual bool element_constraint( bool request_jacobian,
                                     libMesh::DiffContext& context );

    //! Boundary contributions to \f$F(u)\f$ which do not have time varying components.
    virtual bool side_constraint( bool request_jacobian,
                                  libMesh::DiffContext& context );

    //! Contributions to \f$F(u)\f$ on SCALAR variables which do not have time varying components.
    virtual bool nonlocal_constraint( bool request_jacobian,
                                      libMesh::DiffContext& context );

    //! Contributions to \f$C(u)\dot{u}\f$
    virtual bool damping_residual( bool request_jacobian,
                                   libMesh::DiffContext& context );

    //! Contributions to \f$M(u)\dot{u}\f$
    virtual bool mass_residual( bool request_jacobian,
                                libMesh::DiffContext& context );

    //! Contributions to \f$M(u)\dot{u}\f$ on SCALAR variables
    virtual bool nonlocal_mass_residual( bool request_jacobian,
                                         libMesh::DiffContext& context );

    //! Query to check if a particular physics has been enabled
    bool has_physics( const std::string physics_name ) const;

    std::shared_ptr<GRINS::Physics> get_physics( const std::string physics_name );

    std::shared_ptr<GRINS::Physics> get_physics( const std::string physics_name ) const;

    virtual void compute_postprocessed_quantity( unsigned int quantity_index,
                                                 const AssemblyContext& context,
                                                 const libMesh::Point& point,
                                                 libMesh::Real& value );

    std::vector<std::shared_ptr<NeumannBCContainer> >& get_neumann_bcs()
    { return _neumann_bcs; }

    const std::vector<std::shared_ptr<NeumannBCContainer> >& get_neumann_bcs() const
    { return _neumann_bcs; }

  private:

    //! Container of pointers to GRINS::Physics classes requested at runtime.
    /*! Set using the attach_physics_list method as construction is taken care
      of by GRINS::PhysicsFactory. */
    PhysicsList _physics_list;

    bool _use_numerical_jacobians_only;

    // A list of names of variables who need their own numerical
    // jacobian deltas
    std::vector<std::string> _numerical_jacobian_h_variables;

    // A list of values for per-variable numerical jacobian deltas
    std::vector<libMesh::Real> _numerical_jacobian_h_values;

    //! Cached for helping build boundary conditions
    /*! We can't make a copy because it will muck up the UFO detection
      amongst other things. So, we keep a raw pointer. We don't own this
      so we *MUST* not delete. */
    const GetPot* _input;

    //! Neumann boundary conditions
    /*! Store each NeumannBCContainer for each set of BoundaryIDs and Variables,
      as specified in the input file. The container knows what BoundaryIDs and
      Variables it applies to. We use std::shared_ptr here because
      std::unique_ptr may still actually be an AutoPtr. */
    std::vector<std::shared_ptr<NeumannBCContainer> > _neumann_bcs;

    //! Constraint application object
    std::unique_ptr<libMesh::System::Constraint> _constraint;

    // Useful typedef to pointer-to-member functions so we can call all
    // residual and caching functions using a single function (_general_residual)
    typedef void (GRINS::Physics::*ResFuncType) (bool, AssemblyContext &);
    typedef void (GRINS::Physics::*CacheFuncType) (AssemblyContext &);

    // Refactored residual evaluation implementation
    bool _general_residual( bool request_jacobian,
                            libMesh::DiffContext & context,
                            ResFuncType resfunc,
                            CacheFuncType cachefunc);

    //! Extract the bcs from neumann_bcs that are active on bc_id and return them in active_neumann_bcs
    void get_active_neumann_bcs( BoundaryID bc_id,
                                 const std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs,
                                 std::vector<std::shared_ptr<NeumannBCContainer> >& active_neumann_bcs );

    //! Applies the subset of _neumann_bcs that are active on the current element side
    bool apply_neumann_bcs( bool request_jacobian,
                            libMesh::DiffContext& context );
  };

  inline
  std::shared_ptr<GRINS::Physics> MultiphysicsSystem::get_physics( const std::string physics_name ) const
  {
    libmesh_assert(_physics_list.find( physics_name ) != _physics_list.end());

    return _physics_list.find(physics_name)->second;
  }

} //End namespace block

#endif // GRINS_MULTIPHYSICS_SYS_H
