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

#include "reacting_flow_cache.h"

namespace GRINS
{

  ReactingFlowCache::ReactingFlowCache( Real T, Real P, std::vector<Real>& mass_fractions )
    : _T(T),
      _P(P),
      _mass_fractions(mass_fractions),
      _chem_props_set(false),
      _thermo_props_set(false),
      _trans_props_set(false),
      _vel_set(false),
      _vel_grad_set(false),
      _p_hydro_set(false),
      _grad_T_set(false),
      _mf_grad_set(false)
  {
    /*! \todo We need to preallocate the storage here. */
    libmesh_assert_greater(T,0.0);
    libmesh_assert_greater(P,0.0);
    
    return;
  }

  ReactingFlowCache::~ReactingFlowCache()
  {
    return;
  }

  void ReactingFlowCache::set_chemistry_props( Real R, Real M, const std::vector<Real>& omega_dot )
  {
    _R = R;
    _M = M;
    _omega_dot = omega_dot;
    _rho = _P/(_R*_T);
    _chem_props_set = true;
    return;
  }

  void ReactingFlowCache::set_thermo_props( Real cp, std::vector<Real>& h )
  {
    _cp = cp;
    _h = h;
    _thermo_props_set = true;
    return;
  }

  void ReactingFlowCache::set_transport_props( Real mu, Real k, std::vector<Real>& D )
  {
    _mu = mu;
    _k = k;
    _D = D;
    _trans_props_set = true;
    return;
  }

  void ReactingFlowCache::set_velocities( Real u, Real v, Real w )
  {
    _U = libMesh::NumberVectorValue( u, v, w );
    _vel_set = true;
    return;
  }

  void ReactingFlowCache::set_velocity_grads( const libMesh::Gradient& grad_u, const libMesh::Gradient& grad_v,
					      const libMesh::Gradient& grad_w )
  {
    _grad_u = grad_u;
    _grad_v = grad_v;
    _grad_w = grad_w;
    _divU = grad_u(0) + grad_v(1) + grad_w(2);
    _vel_grad_set = true;
    return;
  }

  void ReactingFlowCache::set_p_hydro( Real p_hydro )
  {
    _p_hydro = p_hydro;
    _p_hydro_set = true;
    return;
  }

  void ReactingFlowCache::set_temp_grad( const libMesh::Gradient& grad_T )
  {
    _grad_T = grad_T;
    _grad_T_set = true;
    return;
  }

  void ReactingFlowCache:: set_mass_fractions_grad( const std::vector<libMesh::Gradient>& mass_fractions_grad )
  {
    _mass_fractions_grad = mass_fractions_grad;
    _mf_grad_set = true;
  }

} // namespace GRINS
