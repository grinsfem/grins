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

#include "grins/chemical_species.h"

namespace GRINS
{
  ChemicalSpecies::ChemicalSpecies( const std::string &name, const libMesh::Real mol_wght, const libMesh::Real h_form,
				    const libMesh::Real n_tr_dofs, const int charge )
    : _name      (name),
      _mol_wght  (mol_wght),
      _R         (Constants::R_universal/mol_wght),
      _h_form    (h_form),
      _n_tr_dofs (n_tr_dofs),
      _charge    (charge)
  {
    return;
  }

  ChemicalSpecies::~ChemicalSpecies()
  {
    return;
  }

  void ChemicalSpecies::print (std::ostream &os) const
  {
    os << " -----------------------------\n"
       << "| Species " << this->species() << '\n'
       << " -----------------------------\n"
       << std::scientific
       << "  Mol Wgt = " << this->molar_mass() << '\n'
       << "  R       = " << this->gas_constant()     << '\n'
       << "  h0      = " << this->formation_enthalpy() << '\n'
       << "  n_tr    = " << this->n_tr_dofs() << '\n'
       << "  charge  = " << this->charge() << '\n';

    for (unsigned int l=0; l<this->theta_v().size(); l++)
      os << "  theta_v_" << l << " = " << this->theta_v()[l]
	 << ", ndg = " << this->ndg_v()[l] << "\n";
    
    for (unsigned int l=0; l<this->theta_e().size(); l++)
      os << "  theta_e_" << l << " = " << this->theta_e()[l]
	 << ", ndg = " << this->ndg_e()[l] << "\n";
    
    os << '\n';
    
  }

} //namespace GRINS
