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

#ifndef GRINS_SPALART_ALLMARAS_BC_HANDLING_H
#define GRINS_SPALART_ALLMARAS_BC_HANDLING_H

//GRINS
#include "grins/bc_handling_base.h"

namespace GRINS
{
  // Forward declarations
  class TurbulenceFEVariables;

  //! Base class for reading and handling boundary conditions for physics classes
  class SpalartAllmarasBCHandling : public BCHandlingBase
  {
  public:
    
    SpalartAllmarasBCHandling( const std::string& physics_name, const GetPot& input );
    
    virtual ~SpalartAllmarasBCHandling();

    virtual int string_to_int( const std::string& bc_type_in ) const;

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
    
  protected:

    const TurbulenceFEVariables& _turb_vars;

  private:

    SpalartAllmarasBCHandling();

    enum INS_BC_TYPES{NO_SLIP=1, 
		      PRESCRIBED_VISCOSITY, 
		      PARABOLIC_PROFILE,
		      GENERAL_VISCOSITY};
    
  };
}
#endif // GRINS_SPALART_ALLMARAS_BC_HANDLING_H
