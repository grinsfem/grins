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

#ifndef GRINS_CATALYTIC_WALL_HELPER_H
#define GRINS_CATALYTIC_WALL_HELPER_H

// libMesh
#include "libmesh/libmesh.h"

namespace GRINS
{

  class CatalyticWallHelper
  {
  public:
    
    CatalyticWallHelper( const libMesh::Real R_s, const libMesh::Real M_s,
			 const libMesh::Real gamma_s );

    ~CatalyticWallHelper();
    
     //! \f$ \rho_s \gamma \sqrt{ \frac{R_s T}{2\pi M_s} } \f$
    libMesh::Real omega_dot( const libMesh::Real rho_s, const libMesh::Real T ) const;

    libMesh::Real domega_dot_dws(  const libMesh::Real rho_s, const libMesh::Real w_s,
				   const libMesh::Real T, const libMesh::Real R ) const;

    libMesh::Real domega_dot_dT( const libMesh::Real rho_s, const libMesh::Real T ) const;

    void set_gamma( const libMesh::Real gamma );

  protected:

    libMesh::Real _gamma_s;

    //! \f$ \sqrt{ \frac{R_s}{2\pi M_s} } \f$
    const libMesh::Real _C;

  private:

    CatalyticWallHelper();

  };

  /* ------------------------- Inline Functions -------------------------*/

  inline
  libMesh::Real CatalyticWallHelper::omega_dot( const libMesh::Real rho_s, const libMesh::Real T ) const
  {
    return rho_s*_gamma_s*_C*std::sqrt(T);
  }

  inline
  libMesh::Real CatalyticWallHelper::domega_dot_dws( const libMesh::Real rho_s, const libMesh::Real w_s,
						     const libMesh::Real T, const libMesh::Real R ) const
  {
    return (1.0/w_s - rho_s/R)*(this->omega_dot( rho_s, T ));
  }

  inline
  libMesh::Real CatalyticWallHelper::domega_dot_dT( const libMesh::Real rho_s, const libMesh::Real T ) const
  {
    return -0.5/T*(this->omega_dot( rho_s, T ));
  }

  inline
  void CatalyticWallHelper::set_gamma( const libMesh::Real gamma )
  {
    _gamma_s = gamma;
    return;
  }

} // end namespace GRINS

#endif // GRINS_CATALYTIC_WALL_HELPER_H
