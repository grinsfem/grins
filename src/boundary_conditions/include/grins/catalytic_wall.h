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
#ifndef GRINS_CATALYTIC_WALL_H
#define GRINS_CATALYTIC_WALL_H

// GRINS
#include "grins/math_constants.h"
#include "grins/neumann_func_obj.h"

// libMesh forward declarations
namespace libMesh
{
  class FEMContext;
}

namespace GRINS
{
  // GRINS forward declarations
  class CachedValues;
  class ChemicalMixture;

  class CatalyticWall : public NeumannFuncObj
  {
  public:

    CatalyticWall( const ChemicalMixture& chem_mixture,
		   const unsigned int species_index,
		   const VariableIndex T_var,
		   const libMesh::Real gamma );
    ~CatalyticWall();

    virtual libMesh::Real normal_value( const libMesh::FEMContext& context,
					const CachedValues& cache,
					const unsigned int qp );

    virtual libMesh::Real normal_derivative( const libMesh::FEMContext& context, const CachedValues& cache,
					     const unsigned int qp );

    virtual libMesh::Real normal_derivative( const libMesh::FEMContext& context, const CachedValues& cache,
					     const unsigned int qp, 
					     const GRINS::VariableIndex jac_var );

    libMesh::Real omega_dot( libMesh::Real rho_s, libMesh::Real R_s,
			     libMesh::Real T, libMesh::Real M_s ) const;

    libMesh::Real domega_dot_dws( ) const;

    libMesh::Real domega_dot_dT( ) const;

  protected:

    const ChemicalMixture& _chem_mixture;

    unsigned int _species_index;

    VariableIndex _T_var;

    libMesh::Real _gamma;

  private:

    CatalyticWall();

  };

  /* ------------------------- Inline Functions -------------------------*/

  inline
  libMesh::Real CatalyticWall::omega_dot( libMesh::Real rho_s, libMesh::Real R_s,
					  libMesh::Real T, libMesh::Real M_s ) const
  {
    return rho_s*_gamma*std::sqrt( (R_s*T)/(GRINS::Constants::two_pi*M_s) );
  }

  inline
  libMesh::Real CatalyticWall::domega_dot_dws( ) const
  {
    libmesh_not_implemented();
    return 0.0;
  }

  inline
  libMesh::Real CatalyticWall::domega_dot_dT( ) const
  {
    libmesh_not_implemented();
    return 0.0;
  }

} // namespace GRINS

#endif // GRINS_CATALYTIC_WALL_H
