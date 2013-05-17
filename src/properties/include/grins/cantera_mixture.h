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

#ifndef GRINS_CANTERA_MIXTURE_H
#define GRINS_CANTERA_MIXTURE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_CANTERA

// libMesh
#include "libmesh/threads.h"

// Cantera
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"

// Boost
#include <boost/scoped_ptr.hpp>

// libMesh forward declarations
class GetPot;

namespace
{
  static libMesh::Threads::spin_mutex cantera_mutex;
}

namespace GRINS
{
  class CanteraMixture
  {
  public:

    CanteraMixture( const GetPot& input );
    ~CanteraMixture();

    Cantera::IdealGasMix& get_chemistry();
    const Cantera::IdealGasMix& get_chemistry() const;

    Cantera::Transport& get_transport();

    libMesh::Real M( unsigned int species ) const;

  protected:

    boost::scoped_ptr<Cantera::IdealGasMix> _cantera_gas;

    boost::scoped_ptr<Cantera::Transport> _cantera_transport;

  private:

    CanteraMixture();

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  Cantera::IdealGasMix& CanteraMixture::get_chemistry()
  {
    return (*_cantera_gas);
  }

  inline
  const Cantera::IdealGasMix& CanteraMixture::get_chemistry() const
  {
    return (*_cantera_gas);
  }

  inline
  Cantera::Transport& CanteraMixture::get_transport()
  {
    return (*_cantera_transport);
  }

  inline
  libMesh::Real CanteraMixture::M( unsigned int species ) const
  {
    // Cantera returns molar mass in kg/kmol
    return _cantera_gas->molarMass(species);
  }

} // end namespace GRINS

#endif // GRINS_HAVE_CANTERA

#endif // GRINS_CANTERA_MIXTURE_H
