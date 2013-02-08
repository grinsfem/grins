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

#ifndef GRINS_CHEMICAL_SPECIES_H
#define GRINS_CHEMICAL_SPECIES_H

// C++ includes
#include <string>

// libMesh
#include "libmesh/libmesh.h"

// GRINS
#include "grins/physical_constants.h"

namespace GRINS
{
  //! Class to encapsulate data for each chemical species
  /*!
   * This class is designed to store information relevant to a chemical species.
   * All the data stored is constant for each species, so we store const for each
   * variable. The idea is that this will be placed inside ChemicalMixture, which
   * will be a singleton. This is stolen from the FIN-S class SpeciesChemistry.
   */
  class ChemicalSpecies
  {
  public:
    
    //! Constrctor
    ChemicalSpecies( const std::string &name, 
		     const libMesh::Real        mol_wght,
		     const libMesh::Real        h_form,
		     const libMesh::Real        n_tr_dofs,
		     const int         charge );

    //! Destructor
    ~ChemicalSpecies();

      
    //! Default constructor.
    /*!
     * This is technically required for any
     * std::map value type (or operator[] breaks, at least).  But,
     * we never actually want to create a SpeciesChemistry
     * implicitly, so we throw an error if this is ever used.
     */      
    ChemicalSpecies() : 
      _name("Err"), 
      _mol_wght(0), 
      _R(0),
      _h_form(0), 
      _n_tr_dofs(0), 
      _charge(0)
    { libmesh_error(); }

    //! Returns a descriptive name for this species.
    inline
    const std::string & species() const 
    { return _name; }
      
    //!Returns the molar mass in (g/mol) or (kg/kmol).
    inline
    libMesh::Real molar_mass() const
    { return _mol_wght; }

    //! Returns the species ideal gas constant in [J/kg-K]
    /*!
     * \f$ R \equiv \frac{\hat{R}}{M} \f$ where
     * \f$ R\f$ is the universal gas constant and
     * \f$ M \f$ is the species molar mass.
     */
    inline
    libMesh::Real gas_constant() const 
    { return _R; }

    //! Returns formation enthalpy in units of [J/kg]
    inline
    libMesh::Real formation_enthalpy() const 
    { return _h_form; }
    
    //! Returns number of translational degrees of freedom
    inline
    libMesh::Real n_tr_dofs() const 
    { return _n_tr_dofs; }
    
    //! Returns electrical charge number
    inline
    int charge() const 
    { return _charge; }

    //! Returns true if the chemical species has vibrational degrees 
    inline
    bool has_vibrational_modes() const 
    { return !_theta_v.empty(); }

    //! Returns true if the chemical species has vibrational degrees of freedom.
    inline
    unsigned int n_vibrational_modes() const
    {
      libmesh_assert (_theta_v.size() == _ndg_v.size());
      return _theta_v.size();
    }

    //!Characteristic vibrational temperature [K].
    inline
    const std::vector<libMesh::Real> & theta_v() const
    { return _theta_v; }

    //! Degeneracies for each vibrational mode.
    inline
    const std::vector<unsigned int> & ndg_v() const 
    { return _ndg_v; }

    //! Characteristic electronic excitation temperatures [K].
    inline
    const std::vector<libMesh::Real> & theta_e() const 
    { return _theta_e; }

    //! Degeneracies for each electronic modes.
    inline
    const std::vector<unsigned int> & ndg_e() const
    { return _ndg_e; }

    //! Formatted print
    /*!
     * Defaults to \p std::cout.
     */
    void print(std::ostream &os = std::cout) const;

    //!Formatted print 
    /*!
     * Allows you to do std::cout << object << std::endl;
     */
    friend std::ostream & operator << (std::ostream &os, const ChemicalSpecies &species)
    {
      species.print(os);
      return os;
    }

  protected:

    //! Name of chemical species
    const std::string _name;

    //! Molecular weight (or molar mass) in units of [g/mol] or [kg/kmol]
    const libMesh::Real _mol_wght;

    //! Gas constant in units of [J/kg-K]
    const libMesh::Real _R;

    //! Formation enthalpy in units of [J/kg]
    const libMesh::Real _h_form;

    //! Number of translational degrees of freedom
    const libMesh::Real _n_tr_dofs;

    //! Electrical charge number
    const int _charge;

    //! Characteristic vibrational temperature in units of [K]
    std::vector<libMesh::Real> _theta_v;

    //! Degeneracies for each vibrational mode
    std::vector<unsigned int> _ndg_v;

    //! Characteristic electronic temperature in units of [K]
    std::vector<libMesh::Real> _theta_e;

    //! Degeneracies for each electronic mode
    std::vector<unsigned int> _ndg_e;
    
  }; // class ChemicalSpecies

} //namespace GRINS

#endif // GRINS_CHEMICAL_SPECIES_H
