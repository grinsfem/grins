//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "average_nusselt_number.h"

namespace GRINS
{
  AverageNusseltNumber::AverageNusseltNumber( const GetPot& input )
    : assemble_qoi_sides(true)
  {
    this->read_input_options(input);

    return;
  }

  AverageNusseltNumber::~AverageNusseltNumber()
  {
    return;
  }

  void AverageNusseltNumber::read_input_options( const GetPot& input )
  {
    this->_k = input( "QoI/NusseltNumber/thermal_conductivity", -1.0 );
    if( this->_k < 0.0 )
      {
	std::cerr << "Error: thermal conductivity for AverageNusseltNumber must be positive." << std::endl
		  << "Found k = " << _k << std::endl;
	libmesh_error();
      }

    return;
  }

  void init( const GetPot& input, const libMesh::FEMSystem& system )
  {
    std::string T_var_name = input();
    this->_T_var = system.variable_number(T_var_name);
    return;
  }

  void side_qoi( DiffContext& context, const QoISet& )
  {
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    for( std::set<libMesh::boundary_id_type>::const_iterator id = _bc_ids.begin();
	 id != _bc_ids.end(); id++ )
      {
	if( c.has_boundary_id( (*id) )
	  {
	    // Element Jacobian * quadrature weights for interior integration
	    const std::vector<Real> &JxW = c.element_fe_var[this->_T_var]->get_JxW();
	    
	    unsigned int n_qpoints = (c.get_element_qrule())->n_points();
	    
	    Number& qoi = c.side_qoi[0];
	    
	    // Loop over quadrature points  
	    
	    for (unsigned int qp = 0; qp != n_qpoints; qp++)
	      {
		// Get co-ordinate locations of the current quadrature point
		const Real xf = xyz[qp](0);
		const Real yf = xyz[qp](1);
		
		// If in the sub-domain omega, add the contribution to the integral R
		if(fabs(xf - 0.875) <= 0.125 && fabs(yf - 0.125) <= 0.125)
		  {
		    // Get the solution value at the quadrature point
		    Number T = c.interior_value(this->_T_var, qp);
		    
		    // Update the elemental increment dR for each qp
		    dQoI_0 += JxW[qp] * T;
		  }
		
	      }
	  }
      }

    return;
  }

} //namespace GRINS
