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
#include "grins/mooney_rivlin.h"

// GRINS
#include "grins/common.h"
#include "grins/materials_parsing.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  MooneyRivlin::MooneyRivlin( const GetPot& input )
    : HyperelasticStrainEnergy<MooneyRivlin>(),
    ParameterUser("MooneyRivlin"),
    _C1(-1),
    _C2(-1)
  {
    // Warning about this constructor being deprecated
    {
      std::string warning = "WARNING: Use of this constructor is DEPRECATED.\n";
      warning += "         Please update to use constructor with input material name.\n";
      grins_warning(warning);
    }

    //Force the user to specify C1 and C2
    if( !input.have_variable("Physics/MooneyRivlin/C1") ||
        !input.have_variable("Physics/MooneyRivlin/C2")    )
      {
        std::cerr << "Error: Must specify both C1 and C2 for Mooney-Rivlin material." << std::endl
                  << "       They can be specified in Physics/MooneyRivlin/C1 and Physics/MooneyRivlin/C2" << std::endl;
        libmesh_error();
      }

    this->set_parameter
      (_C1, input, "Physics/MooneyRivlin/C1", _C1);

    this->set_parameter
      (_C2, input, "Physics/MooneyRivlin/C2", _C2);
    return;
  }

  MooneyRivlin::MooneyRivlin( const GetPot& input, const std::string& material )
    : HyperelasticStrainEnergy<MooneyRivlin>(),
    ParameterUser("MooneyRivlin"),
    _C1(-1),
    _C2(-1)
  {
    MaterialsParsing::duplicate_input_test(input,
                                           "Materials/"+material+"/StressStrainLaw/MooneyRivlin/C1",
                                           "Physics/MooneyRivlin/C1");
    MaterialsParsing::duplicate_input_test(input,
                                           "Materials/"+material+"/StressStrainLaw/MooneyRivlin/C2",
                                           "Physics/MooneyRivlin/C2");

    // Parse the new version
    if( input.have_variable("Materials/"+material+"/StressStrainLaw/MooneyRivlin/C1")  &&
        input.have_variable("Materials/"+material+"/StressStrainLaw/MooneyRivlin/C2") )
      {
        this->set_parameter
          (_C1, input, "Materials/"+material+"/StressStrainLaw/MooneyRivlin/C1", _C1);
        this->set_parameter
          (_C2, input, "Materials/"+material+"/StressStrainLaw/MooneyRivlin/C2", _C2);
      }
    // Parse the old version
    else if( input.have_variable("Physics/MooneyRivlin/C1") &&
             input.have_variable("Physics/MooneyRivlin/C2") )
      {
        MaterialsParsing::dep_input_warning( "Physics/MooneyRivlin/C1",
                                             "StressStrainLaw/MooneyRivlin/C1" );
        MaterialsParsing::dep_input_warning( "Physics/MooneyRivlin/C2",
                                             "StressStrainLaw/MooneyRivlin/C2" );

        this->set_parameter
          (_C1, input, "Physics/MooneyRivlin/C1", _C1);
        this->set_parameter
          (_C2, input, "Physics/MooneyRivlin/C2", _C2);
      }
    else
      {
        libmesh_error_msg("ERROR: Could not find consistent Mooney-Rivlin input!");
      }
  }

  MooneyRivlin::~MooneyRivlin()
  {
    return;
  }

} // end namespace GRINS
