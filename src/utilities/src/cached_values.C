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


#include "grins/cached_values.h"

namespace GRINS
{
  CachedValues::CachedValues()
  {
    return;
  }

  CachedValues::~CachedValues()
  {
    return;
  }

  void CachedValues::add_quantity( unsigned int quantity )
  {
    _cache_list.insert(quantity);
    return;
  }

  void CachedValues::add_quantities( const std::set<unsigned int>& cache_list )
  {
    _cache_list.insert(cache_list.begin(), cache_list.end());
    return;
  }

  void CachedValues::clear()
  {
    _cached_values.clear();
    _cached_gradient_values.clear();
    _cached_vector_values.clear();
    _cached_vector_gradient_values.clear();

    return;
  }

  bool CachedValues::is_active(unsigned int quantity)
  {
    bool value = false;

    if( _cache_list.find(quantity) != _cache_list.end() )
      value = true;

    return value;
  }

  void CachedValues::set_values( unsigned int quantity, std::vector<libMesh::Number>& values )
  {
    _cached_values.insert( std::make_pair( quantity, values ) );
    return;
  }

  void CachedValues::set_gradient_values( unsigned int quantity,
                                          std::vector<libMesh::Gradient>& values )
  {
    // Using insert() breaks here. Not entirely sure why...
    _cached_gradient_values[quantity] = values;
    return;
  }

  void CachedValues::set_vector_gradient_values( unsigned int quantity,
                                                 std::vector<std::vector<libMesh::Gradient> >& values )
  {
    _cached_vector_gradient_values[quantity] = values;
    return;
  }

  void CachedValues::set_vector_values( unsigned int quantity, std::vector<std::vector<libMesh::Number> >& values )
  {
    _cached_vector_values[quantity] =  values;
    return;
  }

  const std::vector<libMesh::Number>& CachedValues::get_cached_values( unsigned int quantity ) const
  {
    libmesh_assert( _cached_values.find(quantity) != _cached_values.end() );
    return _cached_values.find(quantity)->second;
  }

  const std::vector<libMesh::Gradient>& CachedValues::get_cached_gradient_values( unsigned int quantity ) const
  {
    libmesh_assert( _cached_gradient_values.find(quantity) != _cached_gradient_values.end() );
    return _cached_gradient_values.find(quantity)->second;
  }

  const std::vector<std::vector<libMesh::Number> >& CachedValues::get_cached_vector_values( unsigned int quantity ) const
  {
    libmesh_assert( _cached_vector_values.find(quantity) != _cached_vector_values.end() );
    return _cached_vector_values.find(quantity)->second;
  }

  const std::vector<std::vector<libMesh::Gradient> >& CachedValues::get_cached_vector_gradient_values( unsigned int quantity ) const
  {
    libmesh_assert( _cached_vector_gradient_values.find(quantity) != _cached_vector_gradient_values.end() );
    return _cached_vector_gradient_values.find(quantity)->second;
  }

} // namespace GRINS
