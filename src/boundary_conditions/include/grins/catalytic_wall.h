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

#ifndef GRINS_CATALYTIC_WALL_H
#define GRINS_CATALYTIC_WALL_H

// Boost
#include "boost/scoped_ptr.hpp"

// GRINS
#include "grins/math_constants.h"
#include "grins/neumann_func_obj.h"
#include "grins/catalycity_base.h"

namespace GRINS
{
  // GRINS forward declarations
  class CachedValues;
  class AssemblyContext;

  template<typename Chemistry>
  class CatalyticWall : public NeumannFuncObj
  {
  public:

    CatalyticWall( const Chemistry& chem_mixture,
		   const unsigned int species_index,
		   CatalycityBase& gamma );

    ~CatalyticWall();

    void init( const VariableIndex T_var );

    virtual libMesh::Real normal_value( const AssemblyContext& context,
					const CachedValues& cache,
					const unsigned int qp );

    virtual libMesh::Real normal_derivative( const AssemblyContext& context, const CachedValues& cache,
					     const unsigned int qp );

    virtual libMesh::Real normal_derivative( const AssemblyContext& context, const CachedValues& cache,
					     const unsigned int qp, 
					     const GRINS::VariableIndex jac_var );

    //! \f$ \rho_s \gamma \sqrt{ \frac{R_s T}{2\pi M_s} } \f$
    libMesh::Real omega_dot( const libMesh::Real rho_s, const libMesh::Real T ) const;

    libMesh::Real domega_dot_dws(  const libMesh::Real rho_s, const libMesh::Real w_s,
				   const libMesh::Real T, const libMesh::Real R ) const;

    libMesh::Real domega_dot_dT( const libMesh::Real rho_s, const libMesh::Real T ) const;

    void set_catalycity_params( const std::vector<libMesh::Real>& params );

  protected:

    const Chemistry& _chemistry;

    unsigned int _species_index;

    boost::scoped_ptr<CatalycityBase> _gamma_s;

    //! \f$ \sqrt{ \frac{R_s}{2\pi M_s} } \f$
    const libMesh::Real _C;

  private:

    CatalyticWall();

  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename Chemistry>
  inline
  libMesh::Real CatalyticWall<Chemistry>::omega_dot( const libMesh::Real rho_s, const libMesh::Real T ) const
  {
    return rho_s*(*_gamma_s)(T)*_C*std::sqrt(T);
  }

  template<typename Chemistry>
  inline
  libMesh::Real CatalyticWall<Chemistry>::domega_dot_dws( const libMesh::Real rho_s, const libMesh::Real w_s,
                                                          const libMesh::Real T, const libMesh::Real R ) const
  {
    return (1.0/w_s - rho_s/R)*(this->omega_dot( rho_s, T ));
  }

  template<typename Chemistry>
  inline
  libMesh::Real CatalyticWall<Chemistry>::domega_dot_dT( const libMesh::Real rho_s, const libMesh::Real T ) const
  {
    libMesh::Real sqrtT = std::sqrt(T);

    return rho_s*_C*( 0.5/sqrtT*(*_gamma_s)(T) + sqrtT*(*_gamma_s).dT(T) );
  }

  template<typename Chemistry>
  inline
  void CatalyticWall<Chemistry>::set_catalycity_params( const std::vector<libMesh::Real>& params )
  {
    _gamma_s->set_params( params );
    return;
  }

} // namespace GRINS

#endif // GRINS_CATALYTIC_WALL_H
