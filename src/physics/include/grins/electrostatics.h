//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_ELECTROSTATICS_H
#define GRINS_ELECTROSTATICS_H

// GRINS
#include "grins/physics.h"

//libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"

// libMesh forward declarations
namespace libMesh
{
  class FEMSystem;
  class FEMContext;
}

namespace GRINS
{

  //! Physics class for Electrostatics
  /*
    This physics class implements electrostatics using the potential form of the
    equations.
  */
  class Electrostatics : public Physics
  {
  public:

    Electrostatics( const std::string& physics_name, const GetPot& input );
    ~Electrostatics();

    //! Initialization  Electrostatics variables
    /*!
      Add velocity and pressure variables to system.
     */
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Sets velocity variables to be time-evolving
    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    // Context initialization
    virtual void init_context( libMesh::FEMContext& context );

    // residual and jacobian calculations
    // element_*, side_* as *time_derivative, *constraint, *mass_residual

    // Time dependent part(s)
    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context,
					  CachedValues& cache );

    virtual void side_time_derivative( bool compute_jacobian,
				       libMesh::FEMContext& context,
				       CachedValues& cache );

    virtual void compute_element_cache( const libMesh::FEMContext& context, 
					const std::vector<libMesh::Point>& points,
					CachedValues& cache );

  protected:

    //! Physical dimension of problem
    /*! \todo Make this static member of base class? */
    unsigned int _dim;

    //! Index for electric potential
    VariableIndex _V_var;

    //! Name of r-velocity
    std::string _V_var_name;

    //! Element type, read from input
    libMeshEnums::FEFamily _V_FE_family;

    //! Temperature element order, read from input
    libMeshEnums::Order _V_order;

    libMesh::Real _sigma;

  private:

    Electrostatics();

  };

} // namespace GRINS

#endif // GRINS_ELECTROSTATICS_H
