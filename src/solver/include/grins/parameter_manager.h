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


#ifndef PARAMETER_MANAGER_H
#define PARAMETER_MANAGER_H

// GRINS

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/parameter_vector.h"

// C++
#include "boost/tr1/memory.hpp"

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;
  class CompositeQoI;



  //! Simple class to hold and initialize a ParameterVector
  /*! Allows the user to specify parameters by name in the input file.
   */
  class ParameterManager
  {
  public:
    
    ParameterManager() {}
    virtual ~ParameterManager() {}

    virtual void initialize( const GetPot& input, 
                             const std::string & parameters_varname,
                             GRINS::MultiphysicsSystem & system,
                             GRINS::CompositeQoI * qoi);

    /*
     * Ordered list of names of independent parameters to study
     */
    std::vector<std::string> parameter_name_list;

    /*
     * Ordered vector of parameter accessors
     */
    libMesh::ParameterVector parameter_vector;
  };

} // end namespace GRINS
#endif // PARAMETER_MANAGER_H
