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

#ifndef GRINS_QOI_OUTPUT_H
#define GRINS_QOI_OUTPUT_H

#include <string>

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  namespace Parallel
  {
    class Communicator;
  }
}

namespace GRINS
{
  // Forward declarations
  class CompositeQoI;

  //! Encapsulate QoI output flags and functionality
  /*! This class handles both the parsing of the options for triggering
    QoI output and implements the functionality for outputting the QoIs
    (which is really just a wrapper around calling output from QoI classes).
    Currently, the user can enable printing the QoI info to the display
    (std::cout) and to a file by specifing the filename in the corresponding
    input option. */
  class QoIOutput
  {
  public:

    QoIOutput( const GetPot & input );

    ~QoIOutput(){}

    //! Function to query whether any input options set to output qoi
    /*! Returns true if user requested to output QoI in any one of the avaiable
      modes, false otherwise. */
    bool output_qoi_set() const
    { return (_output_to_display || _output_to_file ); }

    //! Output the QoI values for all triggered output modes
    /*! This function assumes that the qoi values have been assembled by
      the System. */
    void output_qois( const CompositeQoI & qois, const libMesh::Parallel::Communicator & comm ) const;

  protected:

    bool _output_to_display;

    bool _output_to_file;

    std::string _file_prefix;

  };

} // end namespace GRINS

#endif // GRINS_QOI_OUTPUT_H
