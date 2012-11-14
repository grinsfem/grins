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
      _grad_T_set(false)
  {
    return;
  }

  ReactingFlowCache::~ReactingFlowCache()
  {
    return;
  }

  void ReactingFlowCache::set_chemistry_props( Real R, Real M )
  {
    _R = R;
    _M = M;
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

  void ReactingFlowCache::set_velocity_grads( libMesh::Gradient& grad_u, libMesh::Gradient& grad_v,
					      libMesh::Gradient& grad_w )
  {
    _grad_u = grad_u;
    _grad_v = grad_v;
    _grad_w = grad_w;
    _divU = grad_u(0) + grad_v(1) + grad_w(2);
    _vel_grad_set = true;
    return;
  }

  void ReactingFlowCache::set_p_hyrdo( Real p_hydro )
  {
    _p_hydro = p_hydro;
    _p_hydro_set = true;
    return;
  }

  void ReactingFlowCache::set_temp_grad( libMesh::Gradient& grad_T )
  {
    _grad_T = grad_T;
    _grad_T_set = true;
    return;
  }
}
