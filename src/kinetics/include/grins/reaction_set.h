//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_REACTION_SET_H
#define GRINS_REACTION_SET_H

// libMesh
#include "libmesh/getpot.h"

// GRINS
#include "grins/chemical_mixture.h"
#include "grins/reaction.h"

namespace GRINS
{

   /**
   * This class encapsulates all the reaction mechanisms considered in a
   * chemical nonequilibrium simulation.
   */
   class ReactionSet
   {

   public:

     //! Constructor.
     ReactionSet( const ChemicalMixture& chem_mixture );

     ~ReactionSet();

     //! \returns the number of species.
     unsigned int n_species() const;
     
     //! \returns the number of reactions.
     unsigned int n_reactions() const;
     
     //! Add a reaction to the system.
     void add_reaction(const Reaction& reaction);

     //! \returns a constant reference to reaction \p r.
     const Reaction& reaction(const unsigned int r) const;

     //! Compute species production/destruction rates per unit volume in \f$ \left(kg/sec/m^3\right)\f$
     void compute_mass_sources ( const libMesh::Real T,
				 const libMesh::Real rho,
				 const libMesh::Real R_mix,
				 const std::vector<libMesh::Real>& mass_fractions,
				 const std::vector<libMesh::Real>& molar_densities,
				 const std::vector<libMesh::Real>& h_RT_minus_s_R,
				 std::vector<libMesh::Real>& mass_sources ) const;

     //! Compute the rates of progress for each reaction
     void compute_reaction_rates( const libMesh::Real T,
				  const libMesh::Real rho,
				  const libMesh::Real R_mix,
				  const std::vector<libMesh::Real>& mass_fractions,
				  const std::vector<libMesh::Real>& molar_densities,
				  const std::vector<libMesh::Real>& h_RT_minus_s_R,
				  std::vector<libMesh::Real>& net_reaction_rates ) const;

     //! Formatted print, by default to \p libMesh::out.
     void print( std::ostream& os = libMesh::out ) const;
     
     //! Formatted print.
     friend std::ostream& operator<<( std::ostream& os, const ReactionSet &rset );

   private:
     
     ReactionSet();
     
     const ChemicalMixture& _chem_mixture;

     std::vector<Reaction> _reactions;

   };
  
  /* ------------------------- Inline Functions -------------------------*/
  inline
  unsigned int ReactionSet::n_species() const
  {
    return _chem_mixture.n_species();
  }

  inline
  unsigned int ReactionSet::n_reactions() const
  {
    return _reactions.size();
  }

  inline
  void ReactionSet::add_reaction(const Reaction& reaction)
  {
    _reactions.push_back(reaction);
    
    // and make sure it is initialized!
    _reactions.back().initialize();

    return;
  }
  
  inline
  const Reaction& ReactionSet::reaction(const unsigned int r) const      
  {
    libmesh_assert_less(r, this->n_reactions());
    return _reactions[r];
  }

} // namespace GRINS

#endif // GRINS_REACTION_SET_H
