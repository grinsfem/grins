//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_SOLID_MECHANICS_BC_HANDLING_H
#define GRINS_SOLID_MECHANICS_BC_HANDLING_H

//GRINS
#include "grins/bc_handling_base.h"
#include "grins/solid_mechanics_variables.h"

namespace GRINS
{

  class SolidMechanicsBCHandling : public BCHandlingBase
  {
  public:

    SolidMechanicsBCHandling( const std::string& physics_name,
                              const GetPot& input );

    virtual ~SolidMechanicsBCHandling();

    virtual int string_to_int( const std::string& bc_type_in ) const;

    virtual void init_bc_data( const libMesh::FEMSystem& system );

    virtual void init_bc_types( const GRINS::BoundaryID bc_id,
				const std::string& bc_id_string,
				const int bc_type,
				const std::string& bc_vars,
				const std::string& bc_value,
				const GetPot& input );

    virtual void user_init_dirichlet_bcs( libMesh::FEMSystem* system,
					  libMesh::DofMap& dof_map,
					  GRINS::BoundaryID bc_id,
					  GRINS::BCType bc_type ) const;

    virtual void user_apply_neumann_bcs( AssemblyContext& context,
					 const GRINS::CachedValues& cache,
					 const bool request_jacobian,
					 const GRINS::BoundaryID bc_id,
					 const GRINS::BCType bc_type ) const;

  protected:

    SolidMechanicsVariables _disp_vars;

  private:

    SolidMechanicsBCHandling();

    enum SOLIDS_BC_TYPES{PINNED=1,
                         CONSTANT_DISPLACEMENT,
                         ROLLER_X,
                         ROLLER_Y,
                         ROLLER_Z,
                         CONSTANT_TRACTION};

  };

} // end namespace GRINS

#endif // GRINS_SOLID_MECHANICS_BC_HANDLING_H
