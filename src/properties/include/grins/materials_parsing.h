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

#ifndef GRINS_MATERIALS_PARSING_H
#define GRINS_MATERIALS_PARSING_H

// C++
#include <string>

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  //! Helper functions for parsing material properties
  /*! There's no state needed for these functions so we put them in
      a namespace instead of an object. */
  namespace MaterialsParsing
  {
    //! Check if Physics/physics section has a material variable
    bool have_material( const GetPot& input, const std::string& physics );

    //! Get the name of the material in the Physics/physics section
    void material_name( const GetPot& input, const std::string& physics,
                        std::string& material );

  } // end namespace MaterialsParsing

  inline
  bool MaterialsParsing::have_material( const GetPot& input, const std::string& physics )
  {
    return input.have_variable("Physics/"+physics+"/material");
  }

  inline
  void MaterialsParsing::material_name( const GetPot& input, const std::string& physics,
                                        std::string& material )
  {
    material = input("Physics/"+physics+"/material", "DIE!");
    return;
  }

} // end namespace GRINS

#endif // GRINS_MATERIALS_PARSING_H
