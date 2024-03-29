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
#include "grins/power_law_catalycity.h"

// C++
#include <cmath>

namespace GRINS
{
  PowerLawCatalycity::PowerLawCatalycity( const libMesh::Real gamma0,
                                          const libMesh::Real Tref,
                                          const libMesh::Real alpha )
    : _gamma0(gamma0),
      _Tref(Tref),
      _alpha(alpha)
  {}

  libMesh::Real PowerLawCatalycity::operator()( const libMesh::Real T ) const
  {
    return _gamma0*std::pow( T/_Tref, _alpha);
  }

  libMesh::Real PowerLawCatalycity::dT( const libMesh::Real T ) const
  {
    return (*this)(T)*_alpha/T;
  }

  void PowerLawCatalycity::get_params( std::vector<libMesh::Real> & params )
  {
    libmesh_assert_equal_to( params.size(), 3 );

    params[0] = _gamma0;

    params[1] = _Tref;

    params[2] = _alpha;
  }

  void PowerLawCatalycity::set_params( const std::vector<libMesh::Real>& params )
  {
    libmesh_assert_equal_to( params.size(), 3 );

    _gamma0 = params[0];

    _Tref = params[1];

    _alpha = params[2];
  }

  CatalycityBase* PowerLawCatalycity::clone() const
  {
    return new PowerLawCatalycity( *this );
  }

  void PowerLawCatalycity::set_parameters(const GetPot & input, const std::string & param_base)
  {
    std::string gamma0_str = param_base+"gamma0";
    this->set_parameter(_gamma0,input,gamma0_str,_gamma0);

    std::string Tref_str = param_base+"Tref";
    this->set_parameter(_Tref,input,Tref_str,_Tref);

    std::string alpha_str = param_base+"alpha";
    this->set_parameter(_alpha,input,alpha_str,_alpha);
  }

} // end namespace GRINS
