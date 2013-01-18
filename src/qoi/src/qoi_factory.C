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

#include "grins/qoi_factory.h"

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

  std::tr1::shared_ptr<QoIBase> QoIFactory::build(const GetPot& input)
  {
    /*! \todo Generalize to multiple QoI case when CompositeQoI is implemented in libMesh */
    std::string qoi_name = input("QoI/enabled_qois", "none" );

    std::tr1::shared_ptr<QoIBase> qoi;
    
    if( qoi_name != "none" )
      {
	this->add_qoi( input, qoi_name, qoi );
	
	/*! \todo Generalize to multiple QoI case when CompositeQoI is implemented in libMesh */
	this->check_qoi_physics_consistency( input, qoi_name );
	
	if( input( "screen-options/echo_qoi", false ) )
	  {
	    /*! \todo Generalize to multiple QoI case when CompositeQoI is implemented in libMesh */
	    this->echo_qoi_list( qoi_name );
	  }
      }

    return qoi;
  }

  void QoIFactory::add_qoi( const GetPot& input, const std::string& qoi_name, std::tr1::shared_ptr<QoIBase>& qoi )
  {
    if( qoi_name == avg_nusselt )
      qoi.reset( new AverageNusseltNumber( input ) );

    else if( qoi_name == vorticity )
      qoi.reset( new Vorticity( input ) );

    else
      {
	 libMesh::err << "Error: Invalid QoI name " << qoi_name << std::endl;
	 libmesh_error();
      }

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

  void QoIFactory::echo_qoi_list( const std::string& qoi_name )
  {
    /*! \todo Generalize to multiple QoI case when CompositeQoI is implemented in libMesh */
    std::cout << "==========================================================" << std::endl
	      << "List of Enabled QoIs:" << std::endl
	      << qoi_name << std::endl
	      <<  "==========================================================" << std::endl;
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
