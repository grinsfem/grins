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

#ifndef GRINS_ANTIOCH_MIXTURE_H
#define GRINS_ANTIOCH_MIXTURE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// libMesh
#include "libmesh/libmesh_common.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"
#include "antioch/cea_mixture.h"
#include "antioch/reaction_set.h"

// Boost
#include <boost/scoped_ptr.hpp>

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  class AntiochMixture
  {
  public:

    AntiochMixture( const GetPot& input );
    ~AntiochMixture();

    libMesh::Real M( unsigned int species ) const;

    libMesh::Real M_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real R( unsigned int species ) const;

    libMesh::Real R_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real X( unsigned int species, libMesh::Real M, libMesh::Real mass_fraction ) const;

    void X( libMesh::Real M, const std::vector<libMesh::Real>& mass_fractions, 
	    std::vector<libMesh::Real>& mole_fractions ) const;

    libMesh::Real molar_density( const unsigned int species,
                                 const libMesh::Real rho,
                                 const libMesh::Real mass_fraction ) const;

    void molar_densities( const libMesh::Real rho,
			  const std::vector<libMesh::Real>& mass_fractions,
			  std::vector<libMesh::Real>& molar_densities ) const;

    unsigned int n_species() const;

    unsigned int species_index( const std::string& species_name ) const;

    std::string species_name( unsigned int species_index ) const;

    const Antioch::ChemicalMixture<libMesh::Real>& chemical_mixture() const;

    const Antioch::ReactionSet<libMesh::Real>& reaction_set() const;

    const Antioch::CEAThermoMixture<libMesh::Real>& cea_mixture() const;

  protected:

    boost::scoped_ptr<Antioch::ChemicalMixture<libMesh::Real> > _antioch_gas;

    boost::scoped_ptr<Antioch::ReactionSet<libMesh::Real> > _reaction_set;

    boost::scoped_ptr<Antioch::CEAThermoMixture<libMesh::Real> > _cea_mixture;

  private:

    AntiochMixture();

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  libMesh::Real AntiochMixture::M( unsigned int species ) const 
  {
    return _antioch_gas->M(species);
  }

  inline
  libMesh::Real AntiochMixture::M_mix( const std::vector<libMesh::Real>& mass_fractions ) const 
  {
    return _antioch_gas->M(mass_fractions);
  }

  inline
  libMesh::Real AntiochMixture::R( unsigned int species ) const 
  {
    return _antioch_gas->R(species);
  }

  inline
  libMesh::Real AntiochMixture::R_mix( const std::vector<libMesh::Real>& mass_fractions ) const 
  {
    return _antioch_gas->R(mass_fractions);
  }

  inline
  libMesh::Real AntiochMixture::X( unsigned int species, const libMesh::Real M,
                                   const libMesh::Real mass_fraction ) const
  {
    return _antioch_gas->X(species,M,mass_fraction);
  }

  inline
  void AntiochMixture::X( libMesh::Real M,
                          const std::vector<libMesh::Real>& mass_fractions, 
                          std::vector<libMesh::Real>& mole_fractions ) const
  {
    _antioch_gas->X(M,mass_fractions,mole_fractions);
    return;
  }

  inline
  unsigned int AntiochMixture::n_species() const
  {
    return _antioch_gas->n_species();
  }

  inline
  unsigned int AntiochMixture::species_index( const std::string& species_name ) const
  {
    return _antioch_gas->active_species_name_map().find(species_name)->second;
  }

  inline
  std::string AntiochMixture::species_name( unsigned int /*species_index*/ ) const
  {
    libmesh_not_implemented();
    return "dummy";
  }

  inline
  const Antioch::ChemicalMixture<libMesh::Real>& AntiochMixture::chemical_mixture() const
  {
    return *_antioch_gas.get();
  }

  inline
  const Antioch::ReactionSet<libMesh::Real>& AntiochMixture::reaction_set() const
  {
    return *_reaction_set.get();
  }

  inline
  const Antioch::CEAThermoMixture<libMesh::Real>& AntiochMixture::cea_mixture() const
  {
    return *_cea_mixture.get();
  }

  inline
  libMesh::Real AntiochMixture::molar_density( const unsigned int species,
                                               const libMesh::Real rho,
                                               const libMesh::Real mass_fraction ) const
  {
    return _antioch_gas->molar_density( species, rho, mass_fraction );
  }

  inline
  void AntiochMixture::molar_densities( const libMesh::Real rho,
                                        const std::vector<libMesh::Real>& mass_fractions,
                                        std::vector<libMesh::Real>& molar_densities ) const
  {
    _antioch_gas->molar_densities( rho, mass_fractions, molar_densities );
    return;
  }
  
} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_MIXTURE_H
