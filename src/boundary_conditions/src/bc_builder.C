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
#include "grins/bc_builder.h"

// GRINS
#include "grins/boundary_condition_names.h"
#include "grins/dirichlet_bc_factory_abstract.h"
#include "grins/neumann_bc_factory_abstract.h"
#include "grins/default_bc_builder.h"
#include "grins/old_style_bc_builder.h"

// libMesh
#include "libmesh/dof_map.h"
#include "libmesh/periodic_boundary.h"

namespace GRINS
{
  void BCBuilder::build_boundary_conditions( const GetPot& input,
                                             MultiphysicsSystem& system,
                                             std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs )
  {
    std::unique_ptr<BCBuilder>
      bc_builder = BCBuilder::build_builder(input);

    bc_builder->build_bcs(input,system,neumann_bcs);
  }

  std::unique_ptr<BCBuilder> BCBuilder::build_builder( const GetPot& input )
  {
    bool is_new_style = BCBuilder::is_new_bc_input_style( input );

    std::unique_ptr<BCBuilder> bc_builder;

    if( is_new_style )
      bc_builder.reset( new DefaultBCBuilder );

    else
      bc_builder.reset( new OldStyleBCBuilder );

    return bc_builder;
  }

  bool BCBuilder::is_new_bc_input_style( const GetPot& input )
  {
    bool new_style = false;

    if( input.have_section(BoundaryConditionNames::bc_section()) )
      new_style = true;

    return new_style;
  }

  void BCBuilder::construct_dbc_core( const GetPot& input,
                                      MultiphysicsSystem& system,
                                      const std::set<BoundaryID>& bc_ids,
                                      const FEVariablesBase& fe_var,
                                      const std::string& section,
                                      const std::string& bc_type,
                                      libMesh::DofMap& dof_map )
  {
    // Give the BC factory access to the System
    DirichletBCFactoryAbstract::set_system( system );

    // Give the BC factory access to the GetPot object
    DirichletBCFactoryAbstract::set_getpot(input);

    // Set the boundary id. This gets reset each time inside the factory.
    DirichletBCFactoryAbstract::set_bc_ids( bc_ids );

    // Set the FEVariable. This gets reset each time inside the factory.
    DirichletBCFactoryAbstract::set_fe_var( fe_var );

    // Tell the DirichletBCFactory where to parse the value of the BC,
    // if it needs to
    DirichletBCFactoryAbstract::set_section( section );

    std::unique_ptr<libMesh::DirichletBoundary>
      dbc = DirichletBCFactoryAbstract::build( bc_type );

    dof_map.add_dirichlet_boundary( *dbc );
  }

  void BCBuilder::construct_nbc_core( const GetPot& input,
                                      MultiphysicsSystem& system,
                                      const std::set<BoundaryID>& bc_ids,
                                      const FEVariablesBase& fe_var,
                                      const std::string& section,
                                      const std::string& bc_type,
                                      std::vector<std::shared_ptr<NeumannBCContainer> >& neumann_bcs )
  {
    // Give the BC factory access to the System
    NeumannBCFactoryAbstract::set_system( system );

    // Give the BC factory access to the GetPot object
    NeumannBCFactoryAbstract::set_getpot(input);

    // Set the boundary id. This gets reset each time inside the factory.
    NeumannBCFactoryAbstract::set_bc_ids( bc_ids );

    // Set the FEVariable. This gets reset each time inside the factory.
    NeumannBCFactoryAbstract::set_fe_var( fe_var );

    // Tell the NeumannBCFactory where to parse the value of the BC,
    // if it needs to
    NeumannBCFactoryAbstract::set_section( section );

    std::unique_ptr<NeumannBCContainer>
      nbc = NeumannBCFactoryAbstract::build(bc_type);

    // Get nothing if it's a homogeneous Neumann BC
    // so only try to add it if something was built.
    if(nbc)
      neumann_bcs.push_back(std::move(nbc));
  }

  bool BCBuilder::is_dirichlet_bc_type( const std::string& bc_type )
  {
    // We query the factory to see if it's registered there or not.
    return DirichletBCFactoryAbstract::have_bc_type( bc_type );
  }

  bool BCBuilder::is_neumann_bc_type( const std::string& bc_type )
  {
    // We query the factory to see if it's registered there or not.
    return NeumannBCFactoryAbstract::have_bc_type( bc_type );
  }

  void BCBuilder::add_periodic_bc_to_dofmap( libMesh::boundary_id_type master_id,
                                             libMesh::boundary_id_type slave_id,
                                             const libMesh::RealVectorValue& offset_vector,
                                             libMesh::DofMap& dof_map )
  {
    libMesh::PeriodicBoundary bc( offset_vector );
    bc.myboundary = master_id;
    bc.pairedboundary = slave_id;

    dof_map.add_periodic_boundary( bc );
  }

} // end namespace GRINS
