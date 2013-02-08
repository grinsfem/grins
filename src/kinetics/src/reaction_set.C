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

#include "grins/reaction_set.h"

namespace GRINS
{
  ReactionSet::ReactionSet( const ChemicalMixture& chem_mixture )
    : _chem_mixture(chem_mixture)
  {
    return;
  }

  ReactionSet::~ReactionSet()
  {
    return;
  }

  

  void ReactionSet::compute_mass_sources( const libMesh::Real T,
					  const libMesh::Real rho,
					  const libMesh::Real R_mix,
					  const std::vector<libMesh::Real>& mass_fractions,
					  const std::vector<libMesh::Real>& molar_densities,
					  const std::vector<libMesh::Real>& h_RT_minus_s_R,
					  std::vector<libMesh::Real>& mass_sources ) const
  {
    libmesh_assert_greater(T, 0.0);
    libmesh_assert_greater(rho, 0.0);
    libmesh_assert_greater(R_mix, 0.0);
    libmesh_assert_equal_to( mass_fractions.size(), this->n_species() );
    libmesh_assert_equal_to( molar_densities.size(), this->n_species() );
    libmesh_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );
    libmesh_assert_equal_to( mass_sources.size(), this->n_species() );
    std::fill( mass_sources.begin(), mass_sources.end(), 0.0 );

    // compute the requisite reaction rates
    std::vector<libMesh::Real> net_reaction_rates( this->n_reactions(), 0.0 );
    this->compute_reaction_rates( T, rho, R_mix, mass_fractions, molar_densities,
				  h_RT_minus_s_R, net_reaction_rates );

    // compute the actual mass sources in kmol/sec/m^3
    for (unsigned int rxn = 0; rxn < this->n_reactions(); rxn++)
      {
	const Reaction& reaction = this->reaction(rxn);
	const libMesh::Real rate = net_reaction_rates[rxn];
	
	// reactant contributions
	for (unsigned int r = 0; r  <reaction.n_reactants(); r++)
	  {
	    const unsigned int r_id = reaction.reactant_id(r);
	    const unsigned int r_stoich = reaction.reactant_stoichiometric_coefficient(r);
	    
	    mass_sources[r_id] -= (static_cast<libMesh::Real>(r_stoich)*rate);
	  }
	
	// product contributions
	for (unsigned int p=0; p<reaction.n_products(); p++)
	  {
	    const unsigned int p_id = reaction.product_id(p);
	    const unsigned int p_stoich = reaction.product_stoichiometric_coefficient(p);
	    
	    mass_sources[p_id] += (static_cast<libMesh::Real>(p_stoich)*rate);
	  }
      }

    // finally scale by molar mass
    for (unsigned int s=0; s < this->n_species(); s++)
      {
	mass_sources[s] *= _chem_mixture.M(s);
      }
    
    return;
  }
  
  void ReactionSet::compute_reaction_rates ( const libMesh::Real T,
					     const libMesh::Real rho,
					     const libMesh::Real R_mix,
					     const std::vector<libMesh::Real>& mass_fractions,
					     const std::vector<libMesh::Real>& molar_densities,
					     const std::vector<libMesh::Real>& h_RT_minus_s_R,
					     std::vector<libMesh::Real>& net_reaction_rates ) const
  {
    libmesh_assert_equal_to( net_reaction_rates.size(), this->n_reactions() );
    libmesh_assert_greater(T, 0.0);
    libmesh_assert_greater(rho, 0.0);
    libmesh_assert_greater(R_mix, 0.0);
    libmesh_assert_equal_to( mass_fractions.size(), this->n_species() );
    libmesh_assert_equal_to( molar_densities.size(), this->n_species() );
    libmesh_assert_equal_to( h_RT_minus_s_R.size(), this->n_species() );

    // useful constants
    const libMesh::Real P0    = 1.e5; // standard pressure
    const libMesh::Real RT    = R_mix*T;
    const libMesh::Real P0_RT = P0 / RT; // used to transform equilibrium constant from pressure units

    // compute reaction forward rates & other reaction-sized arrays
    for (unsigned int rxn=0; rxn<this->n_reactions(); rxn++)
      {
	const Reaction& reaction = this->reaction(rxn);

	libMesh::Real kfwd = (reaction.forward_rate())(T);

	libMesh::Real keq = reaction.equilibrium_constant( P0_RT, h_RT_minus_s_R );

	const libMesh::Real kbkwd = kfwd/keq;

	net_reaction_rates[rxn] = reaction.compute_rate_of_progress( molar_densities, kfwd, kbkwd );
      }
    
    return;
  }


  void ReactionSet::print(std::ostream& os) const
  {
    os << "# Number of reactions: " << this->n_reactions() << "\n";

    for (unsigned int r=0; r < this->n_reactions(); r++)
      {
	os << "# " << r << '\n'
	   << this->reaction(r) << "\n";
      }

    return;
  }

  std::ostream& operator<<( std::ostream& os, const ReactionSet &rset )
  {
    rset.print(os);
    return os;
  }
  
} // namespace GRINS
