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
#include "grins/qoi_factory.h"

// GRINS
#include "grins/string_utils.h"
#include "grins/grins_physics_names.h"
#include "grins/qoi_names.h"
#include "grins/average_nusselt_number.h"
#include "grins/vorticity.h"
#include "grins/parsed_boundary_qoi.h"
#include "grins/parsed_interior_qoi.h"

namespace GRINS
{
  QoIFactory::QoIFactory()
  {
    return;
  }
  
  QoIFactory::~QoIFactory()
  {
    return;
  }

  std::tr1::shared_ptr<CompositeQoI> QoIFactory::build(const GetPot& input)
  {
    std::string qoi_list = input("QoI/enabled_qois", "none" );

    std::vector<std::string> qoi_names;

    if( qoi_list != std::string("none") )
      {
        StringUtilities::split_string( qoi_list, std::string(" "), qoi_names );
      }

    std::tr1::shared_ptr<CompositeQoI> qois( new CompositeQoI );
    
    if( !qoi_names.empty() )
      {
        for( std::vector<std::string>::const_iterator name = qoi_names.begin();
             name != qoi_names.end(); ++name )
          {
            this->add_qoi( input, *name, qois );

            this->check_qoi_physics_consistency( input, *name );
          }

	if( input( "screen-options/echo_qoi", false ) )
	  {
	    this->echo_qoi_list( qois );
	  }
      }

    return qois;
  }

  void QoIFactory::add_qoi( const GetPot& /*input*/, const std::string& qoi_name, std::tr1::shared_ptr<CompositeQoI>& qois )
  {
    QoIBase* qoi = NULL;

    if( qoi_name == avg_nusselt )
      {
        qoi = new AverageNusseltNumber( avg_nusselt );
      }

    else if( qoi_name == parsed_boundary )
      {
        qoi =  new ParsedBoundaryQoI( parsed_boundary );
      }

    else if( qoi_name == parsed_interior )
      {
        qoi =  new ParsedInteriorQoI( parsed_interior );
      }

    else if( qoi_name == vorticity )
      {
        qoi =  new Vorticity( vorticity );
      }

    else
      {
	 libMesh::err << "Error: Invalid QoI name " << qoi_name << std::endl;
	 libmesh_error();
      }

    libmesh_assert(qoi);

    qois->add_qoi( *qoi );

    return;
  }

  void QoIFactory::check_qoi_physics_consistency( const GetPot& input, 
						  const std::string& qoi_name )
  {
    int num_physics =  input.vector_variable_size("Physics/enabled_physics");

    // This should be checked other places, but let's be double sure.
    libmesh_assert(num_physics > 0);
  
    std::set<std::string> requested_physics;
    std::set<std::string> required_physics;

    // Build Physics name set
    for( int i = 0; i < num_physics; i++ )
      {
	requested_physics.insert( input("Physics/enabled_physics", "NULL", i ) );
      }
  
    /* If it's Nusselt, we'd better have HeatTransfer or LowMachNavierStokes. 
       HeatTransfer implicitly requires fluids, so no need to check for those. `*/
    if( qoi_name == avg_nusselt )
      {
	required_physics.insert(heat_transfer);
	required_physics.insert(low_mach_navier_stokes);
	this->consistency_helper( requested_physics, required_physics, qoi_name );
      }
      
    return;
  }

  void QoIFactory::echo_qoi_list( std::tr1::shared_ptr<CompositeQoI>& qois )
  {
    /*! \todo Generalize to multiple QoI case when CompositeQoI is implemented in libMesh */
    std::cout << "==========================================================" << std::endl
	      << "List of Enabled QoIs:" << std::endl;
    
    for( unsigned int q = 0; q < qois->n_qois(); q++ )
      {
        std::cout << qois->get_qoi(q).name() << std::endl;      
      }

    std::cout <<  "==========================================================" << std::endl;

    return;
  }

  void QoIFactory::consistency_helper( const std::set<std::string>& requested_physics,
				       const std::set<std::string>& required_physics, 
				       const std::string& qoi_name )
  {
    bool physics_found = false;
    for( std::set<std::string>::const_iterator name = required_physics.begin();
	 name != required_physics.end();
	 name++ )
      {
	if( requested_physics.find( (*name) ) != requested_physics.end() )
	  physics_found = true;
      }

    if( !physics_found )
      this->consistency_error_msg( qoi_name, required_physics );

    return;
  }

  void QoIFactory::consistency_error_msg( const std::string& qoi_name, 
					  const std::set<std::string>& required_physics )
  {
    libMesh::err << "================================================================" << std::endl
		 << "QoI " << qoi_name << std::endl
		 << "requires one of the following physics which were not found:" <<std::endl;
    
    for( std::set<std::string>::const_iterator name = required_physics.begin();
	 name != required_physics.end();
	 name++ )
      {
	libMesh::err << *name << std::endl;
      }
  
    libMesh::err << "================================================================" << std::endl;

    libmesh_error();
  }

} //namespace GRINS
