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

#ifndef GRINS_REACTION_H
#define GRINS_REACTION_H

//C++
#include <string>
#include <vector>

// GRINS
#include "grins/reaction_enum.h"
#include "grins/arrhenius_rate.h"

namespace GRINS
{
  //!A single reaction mechanism. 
  /*!
    This class encapsulates a single reaction mechanism.  The mechanism could be 
    an elementary reaction, or a three-body reaction.  All reactions are assumed
    to be reversible. This class was originally taken from \p FIN-S.
    \todo{Do we want to template this class around the rate type?}
   */
  class Reaction
  {
  public:

    //! Construct a single reaction mechanism.
    Reaction( const unsigned int n_species, const std::string &equation );
    
    ~Reaction();

    unsigned int n_species() const;
    
    //! \returns the equation for this reaction.
    std::string equation() const;
    
    //! Type of reaction. 
    /*! Type of reaction. Presently only ELEMENTARY or THREE_BODY
     *  reversible reactions are considered.
     */
    ReactionType::ReactionType type() const;

    //! Set the type of reaction. 
    /*! Set the type of reaction. Presently only ELEMENTARY or THREE_BODY
     * reversible reactions are considered.
     */
    void set_type( const ReactionType::ReactionType type);

    bool initialized() const;

    //! \returns the number of reactants.
    unsigned int n_reactants() const;

    //! \returns the number of products.
    unsigned int n_products() const;

    //! \returns the name of the \p r th reactant.
    const std::string& reactant_name(const unsigned int r) const;

    //! \returns the name of the \p p th product.
    const std::string& product_name(const unsigned int p) const;

    //!
    unsigned int reactant_id(const unsigned int r) const;

    //!
    unsigned int product_id(const unsigned int p) const;

    //!
    unsigned int reactant_stoichiometric_coefficient(const unsigned int r) const;

    //!
    unsigned int product_stoichiometric_coefficient(const unsigned int p) const;

    //!
    void add_reactant( const std::string &name,
		       const unsigned int r_id,
		       const unsigned int stoichiometric_coeff);

    //!
    void add_product( const std::string &name,
		      const unsigned int p_id,
		      const unsigned int stoichiometric_coeff);

    //!
    void set_efficiency( const std::string &,
			 const unsigned int s,
			 const libMesh::Real efficiency);

    //!
    libMesh::Real efficiency( const unsigned int s) const;

    //! Computes derived quantities.
    void initialize();    

    //!
    int gamma() const;
    
    //!
    libMesh::Real equilibrium_constant( const libMesh::Real P0_RT,
			       const std::vector<libMesh::Real>& h_RT_minus_s_R ) const;

    //!
    void equilibrium_constant_and_derivative( const libMesh::Real T,
					      const libMesh::Real P0_RT,
					      const std::vector<libMesh::Real>& h_RT_minus_s_R,
					      const std::vector<libMesh::Real>& ddT_h_RT_minus_s_R,
					      libMesh::Real& keq,
					      libMesh::Real& dkeq_dT) const;

    //!
    libMesh::Real compute_rate_of_progress( const std::vector<libMesh::Real>& molar_densities,
				   const libMesh::Real kfwd, 
				   const libMesh::Real kbkwd ) const;
    
    //!
    void compute_rate_of_progress_and_derivatives( const std::vector<libMesh::Real>& molar_densities,
						   const std::vector<libMesh::Real>& molar_mass,
						   const libMesh::Real kfwd,  
						   const libMesh::Real dkfwd_dT, 
						   const libMesh::Real kbkwd,
						   const libMesh::Real dkbkwd_dT,
						   libMesh::Real& Rfwd,
						   libMesh::Real& dRfwd_dT,
						   std::vector<libMesh::Real>& dRfwd_drho, 
						   libMesh::Real& Rbkwd,
						   libMesh::Real& dRbkwd_dT,
						   std::vector<libMesh::Real>& dRbkwd_drho) const;

    //! Return const reference to the forward rate object
    const ArrheniusRate& forward_rate() const;

    //! Return writeable reference to the forward rate object
    ArrheniusRate& forward_rate();

    //! Formatted print, by default to \p libMesh::out.
    void print(std::ostream& os = libMesh::out) const;

    //! Formatted print.
    friend std::ostream& operator << (std::ostream& os, const Reaction &rxn);

  private:
    
    unsigned int _n_species;
    ReactionType::ReactionType _type;
    std::string _equation;
    std::vector<std::string> _reactant_names;
    std::vector<std::string> _product_names;
    std::vector<unsigned int> _reactant_ids;
    std::vector<unsigned int> _product_ids;
    std::vector<unsigned int> _reactant_stoichiometry;
    std::vector<unsigned int> _product_stoichiometry;
    std::vector<libMesh::Real>         _efficiencies;
    std::vector<unsigned int> _species_reactant_stoichiometry;
    std::vector<unsigned int> _species_product_stoichiometry;
    std::vector<int>          _species_delta_stoichiometry;
    int _gamma;
    bool _initialized;

    //! The forward reaction rate modified Arrhenius form.
    ArrheniusRate _forward_rate;

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  unsigned int Reaction::n_species() const
  {
    return _n_species;
  }

  inline
  std::string Reaction::equation() const
  {
    return _equation;
  }

  inline
  ReactionType::ReactionType Reaction::type() const
  {
    return _type;
  }

  inline
  void Reaction::set_type( const ReactionType::ReactionType type)
  {
    _type = type;
    return;
  }

  inline
  bool Reaction::initialized() const
  {
    return _initialized;
  }

  inline
  unsigned int Reaction::n_reactants () const
    {
      libmesh_assert_less(_reactant_ids.size(), this->n_species());
      libmesh_assert_equal_to(_reactant_ids.size(), _reactant_stoichiometry.size());
      libmesh_assert_equal_to(_reactant_ids.size(), _reactant_names.size());
      return _reactant_ids.size();
    }

  inline
  unsigned int Reaction::n_products() const
  {
    libmesh_assert_less(_product_ids.size(), this->n_species());
    libmesh_assert_equal_to(_product_ids.size(), _product_stoichiometry.size());
    libmesh_assert_equal_to(_product_ids.size(),  _product_names.size());
    return _product_ids.size();
  }

  inline
  const std::string& Reaction::reactant_name(const unsigned int r) const
  {       
    libmesh_assert_less(r, _reactant_names.size()); 
    return _reactant_names[r]; 
  }

  inline
  const std::string& Reaction::product_name(const unsigned int p) const
  { 
    libmesh_assert_less(p, _product_names.size()); 
    return _product_names[p]; 
  }
  
  inline
  unsigned int Reaction::reactant_id(const unsigned int r) const
  { 
    libmesh_assert_less(r, _reactant_ids.size()); 
    libmesh_assert_less(_reactant_ids[r], this->n_species());
    return _reactant_ids[r]; 
  }

  inline
  unsigned int Reaction::product_id(const unsigned int p) const
  { 
    libmesh_assert_less(p, _product_ids.size()); 
    libmesh_assert_less(_product_ids[p], this->n_species());
    return _product_ids[p]; 
  }

  inline
  unsigned int Reaction::reactant_stoichiometric_coefficient(const unsigned int r) const
  {      
    libmesh_assert_less(r, _reactant_stoichiometry.size());
    libmesh_assert_less(_reactant_ids[r], this->n_species());
    return _reactant_stoichiometry[r];
  }

  inline
  unsigned int Reaction::product_stoichiometric_coefficient(const unsigned int p) const
  {
    libmesh_assert_less(p, _product_stoichiometry.size());
    libmesh_assert_less(_product_ids[p], this->n_species());
    return _product_stoichiometry[p];
  }

  inline
  void Reaction::add_reactant (const std::string &name,
			       const unsigned int r_id,
			       const unsigned int stoichiometric_coeff)
  {
    libmesh_assert_less(r_id, this->n_species());
    _reactant_names.push_back(name);
    _reactant_ids.push_back(r_id);
    _reactant_stoichiometry.push_back(stoichiometric_coeff);
    return;
  }

  inline
  void Reaction::add_product (const std::string &name,
			      const unsigned int p_id,
			      const unsigned int stoichiometric_coeff)
    {
      libmesh_assert_less(p_id, this->n_species());
      _product_names.push_back(name);
      _product_ids.push_back(p_id);
      _product_stoichiometry.push_back(stoichiometric_coeff);
      return;
    }

  inline
  void Reaction::set_efficiency (const std::string &,
				 const unsigned int s,
				 const libMesh::Real efficiency)
  {
    libmesh_assert_less(s, this->n_species());
    libmesh_assert_less(s, _efficiencies.size());
    libmesh_assert_equal_to(_type, ReactionType::THREE_BODY);
    _efficiencies[s] = efficiency;
    return;
  }

  inline
  libMesh::Real Reaction::efficiency( const unsigned int s ) const
  {
    libmesh_assert_less(s, _efficiencies.size());
    libmesh_assert_equal_to(_type, ReactionType::THREE_BODY);
    return _efficiencies[s];
  }

  inline
  int Reaction::gamma () const
  {
    return _gamma;
  }

  inline
  const ArrheniusRate& Reaction::forward_rate() const
  {
    return _forward_rate;
  }

  inline
  ArrheniusRate& Reaction::forward_rate()
  {
    return _forward_rate;
  }
  
} // namespace GRINS

#endif // GRINS_REACTION_H
