//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_CONVECTION_DIFFUSION_H
#define GRINS_CONVECTION_DIFFUSION_H

// GRINS
#include "grins/physics.h"
#include "grins/var_typedefs.h"
#include "grins/grins_enums.h"

//libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/parsed_function.h"

namespace GRINS
{
  class ConvectionDiffusion : public Physics
  {
  public:

    ConvectionDiffusion( const GRINS::PhysicsName& physics_name, const GetPot& input );

    virtual ~ConvectionDiffusion(){};

    //! Initialize variables for this physics.
    virtual void init_variables( libMesh::FEMSystem* system );

    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

    virtual void element_time_derivative( bool compute_jacobian,
					  AssemblyContext& context,
					  CachedValues& cache );

    virtual void mass_residual( bool compute_jacobian,
				AssemblyContext& context,
				CachedValues& cache );

 protected:

    //! Velocity field, \f$ v(x,y,z,t) \f$;
    std::vector<libMesh::ParsedFunction<libMesh::Number> > _v;

    //! Diffusivity, \f$ \kappa(x,y,z,t) \f$
    libMesh::ParsedFunction<libMesh::Number> _kappa;

    std::string _var_name;

    VariableIndex _var;

    //! Element type, read from input
    GRINSEnums::FEFamily _FE_family;

    //! Element orders, read from input
    GRINSEnums::Order _order;

  private:

    ConvectionDiffusion();

  };

} // end namespace GRINS

#endif // GRINS_CONVECTION_DIFFUSION_H
