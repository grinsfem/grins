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


#ifndef GRINS_CACHED_VALUES_H
#define GRINS_CACHED_VALUES_H

//C++
#include <set>
#include <vector>
#include <map>

// libMesh
#include "libmesh/libmesh.h"
#include "libmesh/vector_value.h" //libMesh::Gradient

// GRINS
#include "grins/cached_quantities_enum.h"

namespace GRINS
{
  class CachedValues
  {
  public:

    CachedValues();
    ~CachedValues();

    void add_quantity( unsigned int quantity );

    void add_quantities( const std::set<unsigned int>& cache_list );

    void clear();

    bool is_active(unsigned int quantity);

    void set_values( unsigned int quantity, std::vector<libMesh::Number>& values );

    void set_gradient_values( unsigned int quantity,
			      std::vector<libMesh::Gradient>& values );

    void set_vector_values( unsigned int quantity,
			    std::vector<std::vector<libMesh::Number> >& values );

    void set_vector_gradient_values( unsigned int quantity,
				     std::vector<std::vector<libMesh::Gradient> >& values );

    const std::vector<libMesh::Number>& get_cached_values( unsigned int quantity ) const;
    
    const std::vector<libMesh::Gradient>& get_cached_gradient_values( unsigned int quantity ) const;

    const std::vector<std::vector<libMesh::Number> >& get_cached_vector_values( unsigned int quantity ) const;

    const std::vector<std::vector<libMesh::Gradient> >& get_cached_vector_gradient_values( unsigned int quantity ) const;

  protected:
    
    std::set<unsigned int> _cache_list;

    std::map<unsigned int,std::vector<libMesh::Number> > _cached_values;
    std::map<unsigned int,std::vector<libMesh::Gradient> > _cached_gradient_values;
    std::map<unsigned int,std::vector<std::vector<libMesh::Number> >  > _cached_vector_values;
    std::map<unsigned int,std::vector<std::vector<libMesh::Gradient> >  > _cached_vector_gradient_values;
    
  };

} // namespace GRINS

#endif // GRINS_CACHED_VALUES_H
