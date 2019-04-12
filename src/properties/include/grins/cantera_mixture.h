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

// libMesh
#include "libmesh/auto_ptr.h" // std::unique_ptr

// GRINS
#include "grins/parameter_user.h"

// libMesh forward declarations
class GetPot;

namespace
{
  static libMesh::Threads::spin_mutex cantera_mutex;
}

namespace GRINS
{
  //! Wrapper class for storing state for computing thermochemistry and transport properties using Cantera
  /*!
    This class is expected to be constructed *before* threads have been forked and will
    live during the whole program. Note that this documentation will always
    be built regardless if Cantera is included in the GRINS build or not. Check configure
    output to confirm that Cantera was included in the build if you wish to use it.
  */
  class CanteraMixture : public ParameterUser
  {
  public:

    CanteraMixture( const GetPot& input, const std::string& material );
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

    std::unique_ptr<Cantera::IdealGasMix> _cantera_gas;

    std::unique_ptr<Cantera::Transport> _cantera_transport;

    std::string parse_mixture( const GetPot& input, const std::string& material );

    std::string parse_chem_file( const GetPot& input, const std::string& material );

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
    // Need to return kg/mol to match Antioch
    return _cantera_gas->molecularWeight(species)/1000;
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
