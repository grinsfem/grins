//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "visualization.h"

GRINS::Visualization::Visualization( const GetPot& input )
  : _vis_output_file_prefix( input("vis-options/vis_output_file_prefix", "unknown" ) ),
    _output_residual( input("vis-options/output_residual", false ) )
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

GRINS::Visualization::~Visualization()
{
  return;
}

void GRINS::Solver::output()
{
  this->dump_visualization( this->_vis_output_file_prefix, 0 );

  return;
}

void GRINS::Solver::output( unsigned int time_step )
{
  std::stringstream suffix;

  suffix << time_step;

  std::string filename = this->_vis_output_file_prefix;
  filename+="."+suffix.str();

  this->dump_visualization( filename, time_step );

  return;
}

void GRINS::Solver::dump_visualization( const std::string filename_prefix, const int time_step )
{
  if( this->_vis_output_file_prefix == "unknown" )
    {
      // TODO: Need consisent way to print warning messages.
      std::cout << " WARNING in GRINS::Solver::dump_visualization :" 
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
	  libMesh::TecplotIO(*(this->_mesh),false).write_equation_systems( filename,
									   *(this->_equation_systems) );
	}
      else if ((*format) == "tecplot_binary" ||
	       (*format) == "plt")
	{
	  std::string filename = filename_prefix+".plt";
	  libMesh::TecplotIO(*(this->_mesh),true).write_equation_systems( filename,
									  *(this->_equation_systems) );
	}
      else if ((*format) == "gmv")
	{
	  std::string filename = filename_prefix+".gmv";
	  GMVIO(*(this->_mesh)).write_equation_systems( filename,
						        *(this->_equation_systems) );
	}
      else if ((*format) == "vtu")
	{
	  std::string filename = filename_prefix+".vtu";
	  VTKIO(*(this->_mesh)).write_equation_systems( filename,
						        *(this->_equation_systems) );
	}
      else if ((*format) == "ExodusII")
	{
	  std::string filename = filename_prefix+".exo";
	  
	  // The "1" is hardcoded for the number of time steps because the ExodusII manual states that
	  // it should be the number of timesteps within the file. Here, we are explicitly only doing 
	  // one timestep per file.
	  ExodusII_IO(*(this->_mesh)).write_timestep( filename,
						      *(this->_equation_systems),
						      1,
						      this->_system->time );
	}
      else if ((*format).find("xda") != std::string::npos ||
	       (*format).find("xdr") != std::string::npos)
	{
	  std::string filename = filename_prefix+"."+(*format);
	  const bool binary = ((*format).find("xdr") != std::string::npos);
	  (this->_equation_systems)->write( filename,
					    binary ? libMeshEnums::ENCODE : libMeshEnums::WRITE,
					    EquationSystems::WRITE_DATA | EquationSystems::WRITE_ADDITIONAL_DATA );
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
