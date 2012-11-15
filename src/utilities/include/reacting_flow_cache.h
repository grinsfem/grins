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

#ifndef GRINS_REACTING_FLOW_CACHE_H
#define GRINS_REACTING_FLOW_CACHE_H

// C++
#include <vector>

// libMesh
#include "libmesh_common.h"
#include "vector_value.h"

namespace GRINS
{
  class ReactingFlowCache
  {
  public:

    ReactingFlowCache( Real T, Real P, std::vector<Real>& mass_fractions );
    ~ReactingFlowCache();

    inline
    Real T() const
    { return _T; };

    inline
    const libMesh::Gradient& grad_T() const
    { libmesh_assert(_grad_T_set); return _grad_T; }

    inline
    Real P() const
    { return _P; }

    inline
    Real p_hydro() const
    { libmesh_assert(_p_hydro_set); return _p_hydro;}

    inline
    const std::vector<Real>& mass_fractions() const
    { return _mass_fractions; }

    inline 
    const std::vector<libMesh::Gradient>& mass_fractions_grad() const
    { libmesh_assert(_mf_grad_set); return _mass_fractions_grad; }

    inline
    Real rho() const
    { libmesh_assert(_chem_props_set); return _rho; }

    inline
    Real cp() const
    { libmesh_assert(_thermo_props_set); return _cp; }

    inline
    Real k() const
    { libmesh_assert(_trans_props_set); return _k; }

    inline
    Real mu() const
    { libmesh_assert(_trans_props_set); return _mu; }

    inline
    const std::vector<Real>& D() const
    { libmesh_assert(_trans_props_set); return _D;}

    inline
    const libMesh::NumberVectorValue& U() const
    { libmesh_assert(_vel_set); return _U; }

    inline
    Real divU() const
    { libmesh_assert(_vel_grad_set); return _divU; }

    inline
    const libMesh::Gradient& grad_u() const
    { libmesh_assert(_vel_grad_set); return _grad_u; }

    inline
    const libMesh::Gradient& grad_v() const
    { libmesh_assert(_vel_grad_set); return _grad_v; }

    inline
    const libMesh::Gradient& grad_w() const
    { libmesh_assert(_vel_grad_set); return _grad_w; }

    void set_chemistry_props( Real R, Real M );
    void set_thermo_props( Real cp, std::vector<Real>& h );
    void set_transport_props( Real mu, Real k, std::vector<Real>& D );
    
    void set_velocities( Real u, Real v, Real w = 0.0 );
    void set_velocity_grads( const libMesh::Gradient& grad_u, const libMesh::Gradient& grad_v,
			     const libMesh::Gradient& grad_w );

    void set_p_hydro( Real p_hydro );

    void set_temp_grad( const libMesh::Gradient& grad_T );

    void set_mass_fractions_grad( const std::vector<libMesh::Gradient>& mass_fractions_grad );

  protected:

    const Real _T;
    const Real _P;
    const std::vector<Real> _mass_fractions;
    std::vector<libMesh::Gradient> _mass_fractions_grad;

    Real _rho;

    Real _R;
    Real _M;
    Real _cp;
    std::vector<Real> _h;
    Real _mu;
    Real _k;
    std::vector<Real> _D;

    libMesh::NumberVectorValue _U;
    libMesh::Gradient _grad_u;
    libMesh::Gradient _grad_v;
    libMesh::Gradient _grad_w;
    libMesh::Number _divU;

    Real _p_hydro;

    libMesh::Gradient _grad_T;

    bool _chem_props_set;
    bool _thermo_props_set;
    bool _trans_props_set;
    bool _vel_set;
    bool _vel_grad_set;
    bool _p_hydro_set;
    bool _grad_T_set;
    bool _mf_grad_set;

  private:
    
    ReactingFlowCache();

  };

} // namespace GRINS

#endif //GRINS_REACTING_FLOW_CACHE_H
