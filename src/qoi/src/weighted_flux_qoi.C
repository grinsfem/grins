//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

    // Read boundary ids on which we want to compute fluxes
    int num_bcs =
      input.vector_variable_size("QoI/WeightedFlux/bc_ids");

    if( num_bcs != num_vars )
      {
        std::cerr << "Error: Must specify exactly one boundary id"
                  << " for each specified weighted flux variable."
                  << std::endl
                  << "Found: " << num_bcs << std::endl;
        libmesh_error();
      }

    std::vector<libMesh::boundary_id_type> bc_ids;

    for( int i = 0; i < num_bcs; i++ )
      bc_ids.push_back( input("QoI/WeightedFlux/bc_ids", -1, i ) );

    // Read weight functions with which to compute fluxes
    int num_weights =
      input.vector_variable_size("QoI/WeightedFlux/weights");

    if( num_weights != num_vars )
      {
        std::cerr << "Error: Must specify exactly one weight function"
                  << " for each specified weighted flux variable."
                  << std::endl
                  << "Found: " << num_weights << std::endl;
        libmesh_error();
      }

    for( int i = 0; i < num_weights; i++ )
      {
        const libMesh::boundary_id_type bc_id =
          input("QoI/WeightedFlux/bc_ids", -1, i );

        libmesh_assert_not_equal_to (bc_id, -1);

        const std::string var_name =
          ( input("QoI/WeightedFlux/variables", std::string(""), i ) );

        libmesh_assert_not_equal_to (var_name, std::string(""));

        const std::string func_string =
          input("QoI/WeightedFlux/weights", std::string(""), i);

        libmesh_assert_not_equal_to (func_string, std::string(""));

        std::set<libMesh::boundary_id_type> bc_id_set;
        bc_id_set.insert(bc_id);

        libMesh::ParsedFunction<libMesh::Number> raw_func (func_string);

        std::vector<unsigned int> var_indices;
        var_indices.push_back( system.variable_number( var_name ) );

        libMesh::CompositeFunction<libMesh::Number>* remapped_func =
          new libMesh::CompositeFunction<libMesh::Number>();

        remapped_func->attach_subfunction(raw_func, var_indices);

        libMesh::DirichletBoundary adjoint_bc
          (bc_id_set, var_indices, remapped_func);

        // FIXME: this is an ugly hack
        MultiphysicsSystem & hacked_system =
          const_cast<MultiphysicsSystem&>(system);
        hacked_system.get_dof_map().add_adjoint_dirichlet_boundary
          (adjoint_bc, qoi_num);
      }

    MultiphysicsSystem & hacked_system =
      const_cast<MultiphysicsSystem&>(system);
    hacked_system.reinit_constraints();
  }

} //namespace GRINS
