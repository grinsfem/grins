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


#ifndef GRINS_DISPLACEMENT_FE_VARIABLES_H
#define GRINS_DISPLACEMENT_FE_VARIABLES_H

// GRINS
#include "grins/single_fe_type_variable.h"
#include "grins/displacement_variables.h"

namespace GRINS
{
  class DisplacementFEVariables : public SingleFETypeVariable,
                                  public DisplacementVariables
  {
  public:

    DisplacementFEVariables( const GetPot& input, const std::string& physics_name );
    virtual ~DisplacementFEVariables(){};

    //! Initialize System variables
    /*!
     *  Additional arguments specify whether the spatial mesh is really 2D or 3D.
     * This is needed for cases such as a 1D beam in 2D (is_2D = true) or 3D (is_3D = true)
     * space or 2D shell manifolds in 3D (is_3D = true).
     */
    void init( libMesh::FEMSystem* system, bool is_2D, bool is_3D );

  private:

    DisplacementFEVariables();

  };

} // end namespace GRINS

#endif // GRINS_DISPLACEMENT_FE_VARIABLES_H
