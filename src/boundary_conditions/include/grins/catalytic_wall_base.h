//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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
#include "grins/neumann_bc_abstract.h"
#include <memory>
#include "grins/var_typedefs.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/auto_ptr.h" // std::unique_ptr

namespace libMesh
{
  class FEMSystem;
}

namespace GRINS
{
  // Forward declarations
  class AssemblyContext;
  class CachedValues;

  template<typename Chemistry>
  class CatalyticWallBase : public NeumannBCAbstract
  {
  public:

    CatalyticWallBase( std::shared_ptr<Chemistry>& chem,
                       std::shared_ptr<CatalycityBase>& gamma,
                       const std::vector<VariableIndex>& species_vars,
                       VariableIndex T_var,
                       libMesh::Real p0,
                       unsigned int reactant_species_idx);

    //! Deprecated
    CatalyticWallBase( const Chemistry& chem_mixture,
                       CatalycityBase& gamma,
                       const unsigned int reactant_species_idx );

    virtual ~CatalyticWallBase(){};

    //! Deprecated
    virtual void apply_fluxes( AssemblyContext& context,
                               const CachedValues& cache,
                               const bool request_jacobian ) =0;

    virtual void init( const libMesh::FEMSystem& /*system*/ ){};

    libMesh::Real rho( libMesh::Real T, libMesh::Real p0, libMesh::Real R_mix) const;

    //! \f$ \rho_s \gamma \sqrt{ \frac{R_s T}{2\pi} } \f$
    libMesh::Real omega_dot( const libMesh::Real rho_s, const libMesh::Real T ) const;

    libMesh::Real domega_dot_dws(  const libMesh::Real rho_s, const libMesh::Real w_s,
                                   const libMesh::Real T, const libMesh::Real R ) const;

    libMesh::Real domega_dot_dT( const libMesh::Real rho_s, const libMesh::Real T ) const;

    void set_catalycity_params( const std::vector<libMesh::Real>& params );

    virtual void register_parameter(  const std::string & param_name,
                                      libMesh::ParameterMultiAccessor< libMesh::Number > & param_pointer) const;

  protected:

    //! Temporary helper to deal with intermediate refactoring
    libMesh::Real eval_gamma( libMesh::Real T ) const;

    //! Temporary helper to deal with intermediate refactoring
    libMesh::Real eval_gamma_dT( libMesh::Real T ) const;

    std::shared_ptr<Chemistry> _chem_ptr;

    //! Deprecated
    const Chemistry& _chemistry;

    std::shared_ptr<CatalycityBase> _gamma_ptr;

    //! Deprecated
    std::unique_ptr<CatalycityBase> _gamma_s;

    //! \f$ \sqrt{ \frac{R_s}{2\pi} } \f$
    const libMesh::Real _C;

    //! \todo make const
    std::vector<VariableIndex> _species_vars;

    //! \todo make const
    VariableIndex _T_var;

    //! Thermodynamic pressure
    /*! Currently, we assume that the thermodynamic pressure is constant. This
      is not true in cavity type systems.

      \todo make const */
    libMesh::Real _p0;
  };

  /* ------------------------- Inline Functions -------------------------*/
  template<typename Chemistry>
  inline
  libMesh::Real CatalyticWallBase<Chemistry>::rho( libMesh::Real T, libMesh::Real p0, libMesh::Real R_mix) const
  {
    return p0/(R_mix*T);
  }

  template<typename Chemistry>
  inline
  libMesh::Real CatalyticWallBase<Chemistry>::eval_gamma( libMesh::Real T ) const
  {
    libMesh::Real value;

    if(_gamma_s)
      value = (*_gamma_s)(T);
    else if(_gamma_ptr)
      value = (*_gamma_ptr)(T);
    else
      libmesh_error();

    return value;
  }

  template<typename Chemistry>
  inline
  libMesh::Real CatalyticWallBase<Chemistry>::eval_gamma_dT( libMesh::Real T ) const
  {
    libMesh::Real value;

    if(_gamma_s)
      value = (*_gamma_s).dT(T);
    else if(_gamma_ptr)
      value = (*_gamma_ptr).dT(T);
    else
      libmesh_error();

    return value;
  }

  template<typename Chemistry>
  inline
  libMesh::Real CatalyticWallBase<Chemistry>::omega_dot( const libMesh::Real rho_s, const libMesh::Real T ) const
  {
    return rho_s*this->eval_gamma(T)*_C*std::sqrt(T);
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

    return rho_s*_C*( 0.5/sqrtT*this->eval_gamma(T) + sqrtT*this->eval_gamma_dT(T) );
  }

} // end namespace GRINS

#endif // GRINS_CATALYTIC_WALL_BASE_H
