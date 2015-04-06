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


// This class
#include "grins/weighted_flux_qoi.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/composite_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/parsed_function.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  WeightedFluxQoI::WeightedFluxQoI( const std::string& qoi_name )
    : QoIBase(qoi_name)
  {
    return;
  }

  WeightedFluxQoI::WeightedFluxQoI( const WeightedFluxQoI& original )
    : QoIBase(original.name())
  {
    if (original._adjoint_bc.get())
      this->_adjoint_bc = original._adjoint_bc->clone();
  }


  WeightedFluxQoI::~WeightedFluxQoI()
  {
    return;
  }

  QoIBase* WeightedFluxQoI::clone() const
  {
    return new WeightedFluxQoI( *this );
  }

  void WeightedFluxQoI::init
    (const GetPot& input,
     const MultiphysicsSystem& system,
     unsigned int qoi_num )
  {
    // Read variables for which we want to compute fluxes
    int num_vars = input.vector_variable_size("QoI/WeightedFlux/variables");

    if( num_vars <= 0 )
      {
	std::cerr << "Error: Must specify at least one variable to compute"
		  << " weighted fluxes." << std::endl
		  << "Found: " << num_vars << std::endl;
	libmesh_error();
      }

    for( int i = 0; i < num_vars; i++ )
      {
	_var_names.push_back
          ( input("QoI/WeightedFlux/variables", std::string(""), i ) );
      }

    // Read boundary ids on which we want to compute fluxes
    int num_bcs =  input.vector_variable_size("QoI/WeightedFlux/bc_ids");

    if( num_bcs <= 0 )
      {
	std::cerr << "Error: Must specify at least one boundary id to compute"
		  << " weighted fluxes." << std::endl
		  << "Found: " << num_bcs << std::endl;
	libmesh_error();
      }

    for( int i = 0; i < num_bcs; i++ )
      {
	_bc_ids.insert( input("QoI/WeightedFlux/bc_ids", -1, i ) );
      }


    std::string func_string =
      input( "QoI/WeightedFlux/weights", std::string(""));

    if ( func_string == "" )
      {
	std::cerr << "Error: Must specify a weights function for"
		  << " WeightedFlux QoI" << std::endl;
	libmesh_error();
      }

    libMesh::ParsedFunction<libMesh::Number> raw_func (func_string);

    std::vector<unsigned int> var_indices;
    for( std::vector<VariableName>::const_iterator name = _var_names.begin();
         name != _var_names.end();
         name++ )
      {
        var_indices.push_back( system.variable_number( *name ) );
      }

    libMesh::CompositeFunction<libMesh::Number>* remapped_func =
      new libMesh::CompositeFunction<libMesh::Number>();

    remapped_func->attach_subfunction(raw_func, var_indices);

    _adjoint_bc.reset(remapped_func);

    libMesh::DirichletBoundary adjoint_bc
      (_bc_ids, var_indices, _adjoint_bc.get());

    // FIXME: this is an ugly hack
    MultiphysicsSystem & hacked_system =
      const_cast<MultiphysicsSystem&>(system);
    hacked_system.get_dof_map().add_adjoint_dirichlet_boundary
      (adjoint_bc, qoi_num);
    hacked_system.reinit_constraints();
  }

} //namespace GRINS
