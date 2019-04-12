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
#include "grins/constant_catalycity.h"

namespace GRINS
{
  ConstantCatalycity::ConstantCatalycity( const libMesh::Real gamma )
    : _gamma(gamma)
  {
    return;
  }

  ConstantCatalycity::~ConstantCatalycity()
  {
    return;
  }

  libMesh::Real ConstantCatalycity::operator()( const libMesh::Real /*T*/ ) const
  {
    return _gamma;
  }

  libMesh::Real ConstantCatalycity::dT( const libMesh::Real /*T*/ ) const
  {
    return 0.0;
  }

  void ConstantCatalycity::get_params( std::vector<libMesh::Real> & params )
  {
    libmesh_assert_equal_to( params.size(), 1 );

    params[0] = _gamma;
  }

  void ConstantCatalycity::set_params( const std::vector<libMesh::Real>& params )
  {
    libmesh_assert_equal_to( params.size(), 1 );

    _gamma = params[0];

    return;
  }

  CatalycityBase* ConstantCatalycity::clone() const
  {
    return new ConstantCatalycity( *this );
  }

  void ConstantCatalycity::set_parameters(const GetPot & input, const std::string & param_base)
  {
    std::string gamma_str = param_base+"gamma";
    this->set_parameter(_gamma,input,gamma_str,_gamma);
  }

} // end namespace GRINS
