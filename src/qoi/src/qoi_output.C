//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "grins/qoi_output.h"

// GRINS
#include "grins/common.h"
#include "grins/qoi_options.h"
#include "grins/composite_qoi.h"

// libMesh
#include "libmesh/getpot.h"

// C++
#include <fstream>
#include <iomanip>
#include <sstream>

namespace GRINS
{
  QoIOutput::QoIOutput( const GetPot & input )
    : _output_to_display( input(OutputParsing::output_section()+"/"+OutputParsing::display_section()+"/"+QoIOptions::output_to_display(), false) ),
      _output_to_file( input.have_variable(OutputParsing::output_section()+"/"+QoIOptions::qoi_section()+"/"+QoIOptions::default_file_prefix()) ),
      _file_prefix( input(OutputParsing::output_section()+"/"+QoIOptions::qoi_section()+"/"+QoIOptions::default_file_prefix(), "nofile") )
  {
    if( input.have_variable("screen-options/print_qoi") )
      {
        std::string warning;
        warning = "WARNING: Option screen-options/print_qoi is DEPRECATED!\n";
        warning += "         Please update input to use ";
        warning += OutputParsing::output_section()+"/"+OutputParsing::display_section()+"/"+QoIOptions::output_to_display()+"\n";
        grins_warning(warning);

        _output_to_display = input("screen-options/print_qoi", false);
      }
  }

  void QoIOutput::output_qois( const CompositeQoI & qois, const libMesh::Parallel::Communicator & comm ) const
  {
    if( _output_to_display )
      {
        std::cout << "==========================================================" << std::endl;
        qois.output_qoi( std::cout );
        std::cout << "==========================================================" << std::endl;
      }

    if( _output_to_file )
      {
        if( comm.rank() == 0 )
          {
            std::ofstream output;
            output.open( _file_prefix+".dat" );

            qois.output_qoi(output);

            output.close();
          }
      }
  }
} // end namespace GRINS
