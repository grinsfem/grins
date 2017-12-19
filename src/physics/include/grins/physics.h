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

#ifndef GRINS_PHYSICS_H
#define GRINS_PHYSICS_H

// C++
#include <string>
#include <set>

//GRINS
#include "grins_config.h"
#include "grins/variable_name_defaults.h"
#include "grins/var_typedefs.h"
#include "grins/physics_naming.h"
#include "grins/cached_values.h"
#include "grins/parameter_user.h"
#include "grins/assembly_context.h"

//libMesh
#include "libmesh/libmesh.h"
#include "libmesh/point.h"
#include "libmesh/fe_base.h"
#include "libmesh/system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/auto_ptr.h"

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  template <typename Scalar>
  class CompositeFunction;

  class FEMSystem;
  class Elem;

  template <typename Scalar>
  class ParameterMultiAccessor;
}

//! GRINS namespace
namespace GRINS
{
  // GRINS forward declarations
  class BCHandlingBase;
  class ICHandlingBase;
  class NBCContainer;
  class DBCContainer;
  class MultiphysicsSystem;
  class FEVariablesBase;

  template <typename Scalar>
  class PostProcessedQuantities;

  //! Physics abstract base class. Defines API for physics to be added to MultiphysicsSystem.
  /*!
    This abstract base class defines the API for use within the MultiphysicsSystem. Each physics
    instantiation must initialize its varaibles (and element orders and types), initialize the DiffContext
    for the added variables, and provide element and boundary routines for assembly
    (to be called by MultiphysicsSystem).

    MultiphysicsSystem (through FEMSystem) solves the following equation:

    \f$M(u)\dot{u} = F(u)\f$

    M = mass matrix
    u = solution vector
    F = time derivative

    Note that for the nonlinear system that is solved for implicit
    time stepping is:

    \f$M(u_{\theta})(u^n - u^{n+1}) + \Delta t F(u) = 0\f$
  */

  //TODO: is it F(u) or F(u_{\theta})?

  //  element* routines work on element interiors
  //  side* routines work on element sides

  //  *_time_derivative correspond to calculating terms for F(u)
  //  *_mass_residual correspond to calculating terms for M(u)\dot{u}

  class Physics : public ParameterUser
  {

  public:

    Physics( const GRINS::PhysicsName& physics_name, const GetPot& input );
    virtual ~Physics();

    //! Initialize variables for this physics.
    virtual void init_variables( libMesh::FEMSystem* /*system*/ ){};

    //! Find if current physics is active on supplied element
    virtual bool enabled_on_elem( const libMesh::Elem* elem );

    //! Sets whether this physics is to be solved with a steady solver or not
    /*! Since the member variable is static, only needs to be called on a single
      physics. */
    void set_is_steady( bool is_steady );

    //! Returns whether or not this physics is being solved with a steady solver.
    bool is_steady() const;

    //! Set whether we should treat the problem as axisymmetric
    static void set_is_axisymmetric( bool is_axisymmetric )
    { _is_axisymmetric = is_axisymmetric; }

    static bool is_axisymmetric()
    { return _is_axisymmetric; }

    //! Set which variables are time evolving.
    /*!
      Set those variables which evolve in time (as opposed to variables that behave like constraints).
      This is done separately from init_variables() because the MultiphysicsSystem must initialize
      its base class before time_evolving variables can be set. Default implementation is no
      time evolving variables.
    */
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Any auxillary initialization a Physics class may need
    /*! This is called after all variables are added, so this method can
      safely query the MultiphysicsSystem about variable information. */
    virtual void auxiliary_init( MultiphysicsSystem& system );

    //! Register name of postprocessed quantity with PostProcessedQuantities
    /*!
      Each Physics class will need to cache an unsigned int corresponding to each
      postprocessed quantity. This will be used in computing the values and putting
      them in the CachedVariables object.
    */
    virtual void register_postprocessing_vars( const GetPot& input,
                                               PostProcessedQuantities<libMesh::Real>& postprocessing );

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

    //! Perform any necessary setup before element assembly begins
    virtual void preassembly( MultiphysicsSystem & /*system*/ ){};

    //! Any reinitialization that needs to be done
    /*! This is called through libMesh::FEMSystem::reinit, which is called e.g.
      after adaptive mesh refinement/coarsening. So, for Physics that need to
      reinit internally, then this method should be overridden. */
    virtual void reinit( MultiphysicsSystem & /*system*/ ){};

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool /*compute_jacobian*/,
                                          AssemblyContext & /*context*/ ){}

    //! Time dependent part(s) of physics for boundaries of elements on the domain boundary
    virtual void side_time_derivative( bool /*compute_jacobian*/,
                                       AssemblyContext & /*context*/ ){}

    //! Time dependent part(s) of physics for scalar variables
    virtual void nonlocal_time_derivative( bool /*compute_jacobian*/,
                                           AssemblyContext & /*context*/ ){}

    //! Constraint part(s) of physics for element interiors
    virtual void element_constraint( bool /*compute_jacobian*/,
                                     AssemblyContext & /*context*/ ){}

    //! Constraint part(s) of physics for boundaries of elements on the domain boundary
    virtual void side_constraint( bool /*compute_jacobian*/,
                                  AssemblyContext & /*context*/ ){}

    //! Constraint part(s) of physics for scalar variables
    virtual void nonlocal_constraint( bool /*compute_jacobian*/,
                                      AssemblyContext & /*context*/ ){}

    //! Damping matrix part(s) for element interiors. All boundary terms lie within the time_derivative part
    virtual void damping_residual( bool /*compute_jacobian*/,
                                   AssemblyContext & /*context*/ ){}

    //! Mass matrix part(s) for element interiors. All boundary terms lie within the time_derivative part
    virtual void mass_residual( bool /*compute_jacobian*/,
                                AssemblyContext & /*context*/ ){}

    //! Mass matrix part(s) for scalar variables.
    virtual void nonlocal_mass_residual( bool /*compute_jacobian*/,
                                         AssemblyContext & /*context*/ ){}

    void init_ics( libMesh::FEMSystem* system,
                   libMesh::CompositeFunction<libMesh::Number>& all_ics );

    virtual void compute_element_time_derivative_cache( AssemblyContext & /*context*/ ){}

    virtual void compute_side_time_derivative_cache( AssemblyContext & /*context*/ ){}

    virtual void compute_nonlocal_time_derivative_cache( AssemblyContext & /*context*/ ){}

    virtual void compute_element_constraint_cache( AssemblyContext & /*context*/ ){}

    virtual void compute_side_constraint_cache( AssemblyContext & /*context*/ ){}

    virtual void compute_nonlocal_constraint_cache( AssemblyContext & /*context*/ ){}

    virtual void compute_damping_residual_cache( AssemblyContext & /*context*/ ){}

    virtual void compute_mass_residual_cache( AssemblyContext & /*context*/ ){}

    virtual void compute_nonlocal_mass_residual_cache( AssemblyContext & /*context*/ ){}

    virtual void compute_postprocessed_quantity( unsigned int quantity_index,
                                                 const AssemblyContext& context,
                                                 const libMesh::Point& point,
                                                 libMesh::Real& value );

    ICHandlingBase* get_ic_handler();

  protected:

    /*! \todo This is straight up copied from libMesh. Need to make this available from libMesh. */
    std::unique_ptr<libMesh::FEGenericBase<libMesh::Real> > build_new_fe( const libMesh::Elem* elem,
                                                                          const libMesh::FEGenericBase<libMesh::Real>* fe,
                                                                          const libMesh::Point p );

    void parse_enabled_subdomains( const GetPot& input,
                                   const std::string& physics_name );

    //! Check that var is enabled on at least the subdomains this Physics is
    void check_var_subdomain_consistency( const FEVariablesBase& var ) const;

    //! Name of the physics object. Used for reading physics specific inputs.
    /*! We use a reference because the physics names are const global objects
      in GRINS namespace */
    /*! No, we use a copy, because otherwise as soon as the memory in
     * std::set<std::string> requested_physics gets overwritten we get
     * in trouble. */
    const PhysicsName _physics_name;

    GRINS::ICHandlingBase* _ic_handler;

    //! Subdomains on which the current Physics class is enabled
    std::set<libMesh::subdomain_id_type> _enabled_subdomains;

    //! Caches whether or not the solver that's being used is steady or not.
    /*! This is need, for example, in flow stabilization as the tau terms change
      depending on whether the solver is steady or unsteady. */
    static bool _is_steady;

    //! Caches whether we are solving an axisymmetric problem or not
    static bool _is_axisymmetric;

  private:
    Physics();

  }; // End Physics class declarations

  /* ------------------------- Inline Functions -------------------------*/
  inline
  ICHandlingBase* Physics::get_ic_handler()
  {
    return _ic_handler;
  }

} // End namespace GRINS

#endif //GRINS_PHYSICS_H
