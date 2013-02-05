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

// This class
#include "grins/visualization.h"

// GRINS
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/gmv_io.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/vtk_io.h"

namespace GRINS
{

  Visualization::Visualization( const GetPot& input )
    : _vis_output_file_prefix( input("vis-options/vis_output_file_prefix", "unknown" ) )
  {
    unsigned int num_formats = input.vector_variable_size("vis-options/output_format");

    // If no format specified, default to ExodusII only
    if( num_formats == 0 )
      {
	_output_format.push_back("ExodusII");
      }

    for( unsigned int i = 0; i < num_formats; i++ )
      {
	_output_format.push_back( input("vis-options/output_format", "DIE", i ) );
      }

    return;
  }

  Visualization::~Visualization()
  {
    return;
  }

  void Visualization::output( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system )
  {
    this->dump_visualization( equation_system, _vis_output_file_prefix, 0.0 );

    return;
  }

  void Visualization::output( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
			      const unsigned int time_step, const Real time )
  {
    std::stringstream suffix;

    suffix << time_step;

    std::string filename = this->_vis_output_file_prefix;
    filename+="."+suffix.str();

    this->dump_visualization( equation_system, filename, time );

    return;
  }

  void Visualization::output_residual( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
				       MultiphysicsSystem* system )
  {
    this->output_residual( equation_system, system, 0, 0.0 );
    return;
  }

  void Visualization::dump_visualization( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
					  const std::string& filename_prefix, 
					  const Real time )
  {
    libMesh::MeshBase& mesh = equation_system->get_mesh();

    if( this->_vis_output_file_prefix == "unknown" )
      {
	// TODO: Need consisent way to print warning messages.
	std::cout << " WARNING in Visualization::dump_visualization :" 
		  << " using 'unknown' as file prefix since it was not set " 
		  << std::endl;
      }
  
    for( std::vector<std::string>::const_iterator format = _output_format.begin();
	 format != _output_format.end();
	 format ++ )
      {
	// The following is a modifed copy from the FIN-S code.
	if ((*format) == "tecplot" ||
	    (*format) == "dat")
	  {
	    std::string filename = filename_prefix+".dat";
	    libMesh::TecplotIO(mesh,false).write_equation_systems( filename,
								   *equation_system );
	  }
	else if ((*format) == "tecplot_binary" ||
		 (*format) == "plt")
	  {
	    std::string filename = filename_prefix+".plt";
	    libMesh::TecplotIO(mesh,true).write_equation_systems( filename,
								  *equation_system );
	  }
	else if ((*format) == "gmv")
	  {
	    std::string filename = filename_prefix+".gmv";
	    GMVIO(mesh).write_equation_systems( filename,
						*equation_system );
	  }
	else if ((*format) == "pvtu")
	  {
	    std::string filename = filename_prefix+".pvtu";
	    VTKIO(mesh).write_equation_systems( filename,
						*equation_system );
	  }
	else if ((*format) == "ExodusII")
	  {
	    std::string filename = filename_prefix+".exo";
	  
	    // The "1" is hardcoded for the number of time steps because the ExodusII manual states that
	    // it should be the number of timesteps within the file. Here, we are explicitly only doing 
	    // one timestep per file.
	    ExodusII_IO(mesh).write_timestep( filename,
					      *equation_system,
					      1,
					      time );
	  }
	else if ((*format).find("xda") != std::string::npos ||
		 (*format).find("xdr") != std::string::npos)
	  {
	    std::string filename = filename_prefix+"."+(*format);
	    const bool binary = ((*format).find("xdr") != std::string::npos);
	    equation_system->write( filename,
				    binary ? libMeshEnums::ENCODE : libMeshEnums::WRITE,
				    EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA );
	  }
	else if ((*format) == "mesh_only" )
	  {
	    std::string filename = filename_prefix+"_mesh.xda";
	    equation_system->get_mesh().write( filename );
	  }
	else
	  {
	    // TODO: Do we want to use this to error throughout the code?
	    // TODO: (at least need to pass/print some message/string) - sahni
	    libmesh_error();
	  }
      } // End loop over formats

    return;
  }

} // namespace GRINS
