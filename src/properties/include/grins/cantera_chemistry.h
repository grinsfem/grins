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

#ifndef GRINS_CANTERA_CHEMISTRY_H
#define GRINS_CANTERA_CHEMISTRY_H

#include "grins_config.h"

#ifdef GRINS_HAVE_CANTERA
// Cantera
#include "cantera/IdealGasMix.h"

// Boost
#include <boost/scoped_ptr.hpp>

// libMesh
#include "libmesh/libmesh_common.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  class CanteraChemistry
  {
  public:

    CanteraChemistry( const GetPot& input );
    ~CanteraChemistry();

    libMesh::Real M( unsigned int species ) const;

    libMesh::Real M_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real R( unsigned int species ) const;

    libMesh::Real R_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real X( unsigned int species, libMesh::Real M, libMesh::Real mass_fraction ) const;

    void X( libMesh::Real M, const std::vector<libMesh::Real>& mass_fractions, 
	    std::vector<libMesh::Real>& mole_fractions ) const;

  protected:

    boost::scoped_ptr<Cantera::IdealGasMix> _cantera_gas;

  private:

    CanteraChemistry();
  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  libMesh::Real CanteraChemistry::M( unsigned int species ) const
  {
    // Cantera returns molar mass in kg/kmol
    return _cantera_gas->molarMass(species);
  }

  inline
  libMesh::Real CanteraChemistry::R( unsigned int species ) const
  {
    // Cantera::GasConstant in J/kmol-K
    // Cantera returns molar mass in kg/kmol
    return Cantera::GasConstant/_cantera_gas->molarMass(species);
  }

  inline
  libMesh::Real CanteraChemistry::X( unsigned int species,
                                     libMesh::Real M_mix,
                                     libMesh::Real mass_fraction ) const
  {
    return mass_fraction*M_mix/this->M(species);
  }

} // end namespace GRINS

#endif // GRINS_HAVE_CANTERA

#endif // GRINS_CANTERA_CHEMISTRY_H
