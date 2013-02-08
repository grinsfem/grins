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
#ifndef GRINS_INC_NAVIER_STOKES_BC_HANDLING_H
#define GRINS_INC_NAVIER_STOKES_BC_HANDLING_H

//GRINS
#include "grins/bc_handling_base.h"

namespace GRINS
{
  //! Base class for reading and handling boundary conditions for physics classes
  class IncompressibleNavierStokesBCHandling : public BCHandlingBase
  {
  public:
    
    IncompressibleNavierStokesBCHandling( const std::string& physics_name,
					  const GetPot& input );
    
    virtual ~IncompressibleNavierStokesBCHandling();

    virtual int string_to_int( const std::string& bc_type_in ) const;

    virtual void init_bc_data( const libMesh::FEMSystem& system );
    
    virtual void init_bc_types( const GRINS::BoundaryID bc_id, 
				const std::string& bc_id_string, 
				const int bc_type, 
				const GetPot& input );

    virtual void user_init_dirichlet_bcs( libMesh::FEMSystem* system,
					  libMesh::DofMap& dof_map,
					  GRINS::BoundaryID bc_id,
					  GRINS::BCType bc_type ) const;
    
  protected:

    std::string _u_var_name, _v_var_name, _w_var_name;

    VariableIndex _u_var, _v_var, _w_var;

  private:

    IncompressibleNavierStokesBCHandling();

    enum INS_BC_TYPES{NO_SLIP=1, 
		      PRESCRIBED_VELOCITY, 
		      PARABOLIC_PROFILE,
		      GENERAL_VELOCITY};

  };

} // end namespace GRINS
#endif // GRINS_INC_NAVIER_STOKES_BC_HANDLING_H
