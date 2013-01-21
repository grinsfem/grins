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

#include "grins/cantera_singleton.h"

namespace GRINS
{
#ifdef GRINS_HAVE_CANTERA
  //Note the shared_ptr default constructor is such that it looks like a NULL pointer.
  std::tr1::shared_ptr<Cantera::IdealGasMix> CanteraSingleton::_cantera = 
    std::tr1::shared_ptr<Cantera::IdealGasMix>();

  Cantera::IdealGasMix& CanteraSingleton::cantera_instance( const GetPot& input )
  {
    // Pointer is null, so we create an instance.
    if( !_cantera )
      {
	const std::string cantera_chem_file = input( "Physics/Chemistry/chem_file", "DIE!" );
	const std::string mixture = input( "Physics/Chemistry/mixture", "DIE!" );

	try
	  {
	    _cantera.reset( new Cantera::IdealGasMix( cantera_chem_file, mixture ) );
	  }
	catch(Cantera::CanteraError)
	  {
	    Cantera::showErrors(std::cerr);
	    libmesh_error();
	  }
      }
    
    // Return a reference
    return *(_cantera.get());
  }

  Cantera::IdealGasMix& CanteraSingleton::cantera_instance()
  {
    if( !_cantera )
      {
	std::cerr << "Error: Must first initialize cantera by calling cantera_instance( const GetPot& input ) " << std::endl;
	libmesh_error();
      }

    return *(_cantera.get());
  }

#endif //GRINS_HAVE_CANTERA
}
