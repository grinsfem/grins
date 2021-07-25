//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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

// This class
#include "grins/arrhenius_catalycity.h"

// C++
#include <cmath>

namespace GRINS
{
  ArrheniusCatalycity::ArrheniusCatalycity( const libMesh::Real gamma0,
                                            const libMesh::Real Ta )
    : _gamma0(gamma0),
      _Ta(Ta)
  {
    return;
  }

  libMesh::Real ArrheniusCatalycity::operator()( const libMesh::Real T ) const
  {
    return _gamma0*std::exp(-_Ta/T);
  }

  libMesh::Real ArrheniusCatalycity::dT( const libMesh::Real T ) const
  {
    return _gamma0*_Ta/(T*T)*std::exp(-_Ta/T);
  }

  void ArrheniusCatalycity::get_params( std::vector<libMesh::Real> & params )
  {
    libmesh_assert_equal_to( params.size(), 2 );

    params[0] = _gamma0;

    params[1] = _Ta;
  }

  void ArrheniusCatalycity::set_params( const std::vector<libMesh::Real>& params )
  {
    libmesh_assert_equal_to( params.size(), 2 );

    _gamma0 = params[0];

    _Ta = params[1];

    return;
  }

  CatalycityBase* ArrheniusCatalycity::clone() const
  {
    return new ArrheniusCatalycity( *this );
  }

  void ArrheniusCatalycity::set_parameters(const GetPot & input, const std::string & param_base)
  {
    std::string gamma0_str = param_base+"gamma0";
    this->set_parameter(_gamma0,input,gamma0_str,_gamma0);

    std::string Ta_str = param_base+"Ta";
    this->set_parameter(_Ta,input,Ta_str,_Ta);
  }

} // end namespace GRINS
