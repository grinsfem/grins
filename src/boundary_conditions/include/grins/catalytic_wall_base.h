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

#ifndef GRINS_CATALYTIC_WALL_BASE_H
#define GRINS_CATALYTIC_WALL_BASE_H

// GRINS
#include "grins/catalycity_base.h"

// libMesh
#include "libmesh/libmesh_common.h"
namespace libMesh
{
  class FEMSystem;
}

namespace GRINS
{

  template<typename Chemistry>
  class CatalyticWallBase
  {
  public:

    CatalyticWallBase( const Chemistry& chem_mixture,
                       CatalycityBase& gamma,
                       const unsigned int reactant_species_idx );

    virtual ~CatalyticWallBase();

    virtual void apply_fluxes( AssemblyContext& context,
                               const CachedValues& cache,
                               const bool request_jacobian ) =0;

    virtual void init( const libMesh::FEMSystem& system );

    void set_axisymmetric( bool is_axisymmetric );

    //! \f$ \rho_s \gamma \sqrt{ \frac{R_s T}{2\pi} } \f$
    libMesh::Real omega_dot( const libMesh::Real rho_s, const libMesh::Real T ) const;

    libMesh::Real domega_dot_dws(  const libMesh::Real rho_s, const libMesh::Real w_s,
				   const libMesh::Real T, const libMesh::Real R ) const;

    libMesh::Real domega_dot_dT( const libMesh::Real rho_s, const libMesh::Real T ) const;

    void set_catalycity_params( const std::vector<libMesh::Real>& params );

  protected:

    const Chemistry& _chemistry;

    boost::scoped_ptr<CatalycityBase> _gamma_s;

    //! \f$ \sqrt{ \frac{R_s}{2\pi} } \f$
    const libMesh::Real _C;

    bool _is_axisymmetric;

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename Chemistry>
  inline
  libMesh::Real CatalyticWallBase<Chemistry>::omega_dot( const libMesh::Real rho_s, const libMesh::Real T ) const
  {
    return rho_s*(*_gamma_s)(T)*_C*std::sqrt(T);
  }

  template<typename Chemistry>
  inline
  libMesh::Real CatalyticWallBase<Chemistry>::domega_dot_dws( const libMesh::Real rho_s, const libMesh::Real w_s,
                                                          const libMesh::Real T, const libMesh::Real R ) const
  {
    return (1.0/w_s - rho_s/R)*(this->omega_dot( rho_s, T ));
  }

  template<typename Chemistry>
  inline
  libMesh::Real CatalyticWallBase<Chemistry>::domega_dot_dT( const libMesh::Real rho_s, const libMesh::Real T ) const
  {
    libMesh::Real sqrtT = std::sqrt(T);

    return rho_s*_C*( 0.5/sqrtT*(*_gamma_s)(T) + sqrtT*(*_gamma_s).dT(T) );
  }

} // end namespace GRINS

#endif // GRINS_CATALYTIC_WALL_BASE_H
