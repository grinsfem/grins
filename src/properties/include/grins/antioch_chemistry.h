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


#ifndef GRINS_ANTIOCH_CHEMISTRY_H
#define GRINS_ANTIOCH_CHEMISTRY_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/property_types.h"
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/auto_ptr.h" // libMesh::UniquePtr

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  //! Wrapper class for Antioch::ChemicalMixture
  /*!
    This class is expected to be constructed *before* threads have been forked and will
    live during the whole program.
    By default, Antioch is working in SI units. Note that this documentation will always
    be built regardless if Antioch is included in the GRINS build or not. Check configure
    output to confirm that Antioch was included in the build.
   */
  class AntiochChemistry : public ParameterUser
  {
  public:

    AntiochChemistry( const GetPot& input, const std::string& material );

    virtual ~AntiochChemistry();

    //! Species molar mass (molecular weight), [kg/mol]
    libMesh::Real M( unsigned int species ) const;

    //! Mixture molar mass (molecular weight), [kg/mol]
    libMesh::Real M_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    //! Species gas constant, [J/kg-K]
    /*! R_universal/M(species) */
    libMesh::Real R( unsigned int species ) const;

    //! Mixture gas constant, [J/kg-K]
    libMesh::Real R_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    //! Species mole fraction, unitless
    libMesh::Real X( unsigned int species, libMesh::Real M, libMesh::Real mass_fraction ) const;

    //! Mole fraction for all species, unitless
    void X( libMesh::Real M, const std::vector<libMesh::Real>& mass_fractions,
	    std::vector<libMesh::Real>& mole_fractions ) const;

    //! Species molar density, [mol/m^3]
    libMesh::Real molar_density( const unsigned int species,
                                 const libMesh::Real rho,
                                 const libMesh::Real mass_fraction ) const;

    //! Molar density for all species, [mol/m^3]
    void molar_densities( const libMesh::Real rho,
			  const std::vector<libMesh::Real>& mass_fractions,
			  std::vector<libMesh::Real>& molar_densities ) const;

    unsigned int n_species() const;

    unsigned int species_index( const std::string& species_name ) const;

    std::string species_name( unsigned int species_index ) const;

    //! Accessor to underlying Antioch object
    const Antioch::ChemicalMixture<libMesh::Real>& chemical_mixture() const;

    //! Accessor for this class
    const AntiochChemistry& chemistry() const;

  protected:

    libMesh::UniquePtr<Antioch::ChemicalMixture<libMesh::Real> > _antioch_gas;

  private:

    AntiochChemistry();

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  libMesh::Real AntiochChemistry::M( unsigned int species ) const
  {
    return _antioch_gas->M(species);
  }

  inline
  libMesh::Real AntiochChemistry::M_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    return _antioch_gas->M(mass_fractions);
  }

  inline
  libMesh::Real AntiochChemistry::R( unsigned int species ) const
  {
    return _antioch_gas->R(species);
  }

  inline
  libMesh::Real AntiochChemistry::R_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    return _antioch_gas->R(mass_fractions);
  }

  inline
  libMesh::Real AntiochChemistry::X( unsigned int species, const libMesh::Real M,
                                     const libMesh::Real mass_fraction ) const
  {
    return _antioch_gas->X(species,M,mass_fraction);
  }

  inline
  void AntiochChemistry::X( libMesh::Real M,
                            const std::vector<libMesh::Real>& mass_fractions,
                            std::vector<libMesh::Real>& mole_fractions ) const
  {
    _antioch_gas->X(M,mass_fractions,mole_fractions);
    return;
  }

  inline
  unsigned int AntiochChemistry::n_species() const
  {
    return _antioch_gas->n_species();
  }

  inline
  unsigned int AntiochChemistry::species_index( const std::string& species_name ) const
  {
    return _antioch_gas->species_name_map().find(species_name)->second;
  }

  inline
  const Antioch::ChemicalMixture<libMesh::Real>& AntiochChemistry::chemical_mixture() const
  {
    return *_antioch_gas.get();
  }

  inline
  libMesh::Real AntiochChemistry::molar_density( const unsigned int species,
                                                 const libMesh::Real rho,
                                                 const libMesh::Real mass_fraction ) const
  {
    return _antioch_gas->molar_density( species, rho, mass_fraction );
  }

  inline
  void AntiochChemistry::molar_densities( const libMesh::Real rho,
                                          const std::vector<libMesh::Real>& mass_fractions,
                                          std::vector<libMesh::Real>& molar_densities ) const
  {
    _antioch_gas->molar_densities( rho, mass_fractions, molar_densities );
    return;
  }

  inline
  const AntiochChemistry& AntiochChemistry::chemistry() const
  {
    return *this;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_CHEMISTRY_H
