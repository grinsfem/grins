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

#ifndef GRINS_CACHED_VALUES_H
#define GRINS_CACHED_VALUES_H

//C++
#include <set>
#include <vector>
#include <map>

// libMesh
#include "libmesh.h"

// GRINS
#include "cached_quantities_enum.h"

namespace GRINS
{
  class CachedValues
  {
  public:

    CachedValues();
    ~CachedValues();

    void add_quantity( Cache::CachedQuantities quantity );

    void add_quantities( const std::set<Cache::CachedQuantities>& cache_list );

    void clear();

    bool is_active(Cache::CachedQuantities quantity);

    void set_values( Cache::CachedQuantities quantity, std::vector<Number>& values );

    void set_vector_values( Cache::CachedQuantities quantity, std::vector<std::vector<Number> >& values );

    const std::vector<Number>& get_cached_values( Cache::CachedQuantities quantity ) const;

    const std::vector<std::vector<Number> >& get_cached_vector_values( Cache::CachedQuantities quantity ) const;

  protected:
    
    std::set<Cache::CachedQuantities> _cache_list;

    std::map<Cache::CachedQuantities,std::vector<Number> > _cached_values;
    std::map<Cache::CachedQuantities,std::vector<std::vector<Number> >  > _cached_vector_values;
    
  };

} // namespace GRINS

#endif // GRINS_CACHED_VALUES_H
