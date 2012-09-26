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
    : QoIBase()
  {
    this->assemble_qoi_sides = true;
    this->assemble_qoi_elements = false;
    this->read_input_options(input);

    return;
  }

  AverageNusseltNumber::~AverageNusseltNumber()
  {
    return;
  }

  libMesh::AutoPtr<libMesh::DifferentiableQoI> AverageNusseltNumber::clone()
  {
    return libMesh::AutoPtr<libMesh::DifferentiableQoI>( new AverageNusseltNumber( *this ) );
  }

  void AverageNusseltNumber::read_input_options( const GetPot& input )
  {
    // Read thermal conductivity
    this->_k = input( "QoI/NusseltNumber/thermal_conductivity", -1.0 );

    if( this->_k < 0.0 )
      {
	std::cerr << "Error: thermal conductivity for AverageNusseltNumber must be positive." << std::endl
		  << "Found k = " << _k << std::endl;
	libmesh_error();
      }

    // Read boundary ids for which we want to compute
    int num_bcs =  input.vector_variable_size("QoI/NusseltNumber/bc_ids");

    if( num_bcs <= 0 )
      {
	std::cerr << "Error: Must specify at least one boundary id to compute"
		  << " average Nusselt number." << std::endl
		  << "Found: " << num_bcs << std::endl;
	libmesh_error();
      }

    for( int i = 0; i < num_bcs; i++ )
      {
	_bc_ids.insert( input("QoI/NusseltNumber/bc_ids", -1, i ) );
      }

    this->_scaling = input( "QoI/NusseltNumber/scaling", 1.0 );

    return;
  }

  void AverageNusseltNumber::init( const GetPot& input, const libMesh::FEMSystem& system )
  {
    // Grab temperature variable index
    std::string T_var_name = input("Physics/VariableNames/Temperature",
				   T_var_name_default);

    this->_T_var = system.variable_number(T_var_name);
    return;
  }

  void AverageNusseltNumber::side_qoi( DiffContext& context, const QoISet& )
  {
    FEMContext &c = libmesh_cast_ref<FEMContext&>(context);

    for( std::set<libMesh::boundary_id_type>::const_iterator id = _bc_ids.begin();
	 id != _bc_ids.end(); id++ )
      {
	if( c.has_side_boundary_id( (*id) ) )
	  {
	    FEBase* side_fe;
	    c.get_side_fe<Real>(this->_T_var, side_fe);

	    const std::vector<Real> &JxW = side_fe->get_JxW();
	    
	    const std::vector<Point>& normals = side_fe->get_normals();

	    unsigned int n_qpoints = (c.get_side_qrule())->n_points();
	    
	    Number& qoi = c.elem_qoi[0];
	    
	    // Loop over quadrature points  
	    
	    for (unsigned int qp = 0; qp != n_qpoints; qp++)
	      {
		// Get the solution value at the quadrature point
		Gradient grad_T = 0.0; 
		c.side_gradient<Real>(this->_T_var, qp, grad_T);
		
		// Update the elemental increment dR for each qp
		qoi += (this->_scaling)*(this->_k)*(grad_T*normals[qp])*JxW[qp];
	      } // quadrature loop

	  } // end check on boundary id

      }

    return;
  }

} //namespace GRINS
