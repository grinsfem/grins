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


#ifndef GRINS_CANTERA_MIXTURE_H
#define GRINS_CANTERA_MIXTURE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_CANTERA

// libMesh
#include "libmesh/threads.h"

// Cantera (with compiler warnings disabled)
#include "libmesh/ignore_warnings.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "libmesh/restore_warnings.h"

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

    libMesh::Real M_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real R( unsigned int species ) const;

    libMesh::Real R_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real X( unsigned int species, libMesh::Real M, libMesh::Real mass_fraction ) const;

    void X( libMesh::Real M, const std::vector<libMesh::Real>& mass_fractions, 
	    std::vector<libMesh::Real>& mole_fractions ) const;

    unsigned int n_species() const;

    unsigned int species_index( const std::string& species_name ) const;

    std::string species_name( unsigned int species_index ) const;

    const CanteraMixture& chemistry() const;

    //! This is basically dummy, but is needed for template games elsewhere.
    typedef CanteraMixture ChemistryParent;

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
    return _cantera_gas->molecularWeight(species);
  }

  inline
  libMesh::Real CanteraMixture::R( unsigned int species ) const
  {
    // Cantera::GasConstant in J/kmol-K
    // Cantera returns molar mass in kg/kmol
    return Cantera::GasConstant/_cantera_gas->molecularWeight(species);
  }

  inline
  libMesh::Real CanteraMixture::X( unsigned int species,
                                   libMesh::Real M_mix,
                                   libMesh::Real mass_fraction ) const
  {
    return mass_fraction*M_mix/this->M(species);
  }

  inline
  unsigned int CanteraMixture::n_species() const
  {
    return _cantera_gas->nSpecies();
  }

  inline
  unsigned int CanteraMixture::species_index( const std::string& species_name ) const
  {
    return _cantera_gas->speciesIndex( species_name );
  }

  inline
  std::string CanteraMixture::species_name( unsigned int species_index ) const
  {
    return _cantera_gas->speciesName( species_index );
  }

  inline
  const CanteraMixture& CanteraMixture::chemistry() const
  {
    return *this;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_CANTERA

#endif // GRINS_CANTERA_MIXTURE_H
