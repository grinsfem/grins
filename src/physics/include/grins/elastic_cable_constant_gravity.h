//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_ELASTIC_CABLE_CONSTANT_GRAVITY_H
#define GRINS_ELASTIC_CABLE_CONSTANT_GRAVITY_H

//GRINS
#include "grins/elastic_cable_base.h"

namespace GRINS
{
  class ElasticCableConstantGravity : public ElasticCableBase
  {
  public:
	  ElasticCableConstantGravity( const GRINS::PhysicsName& physics_name, const GetPot& input );
    virtual ~ElasticCableConstantGravity();

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context,
                                          CachedValues& cache );

  private:

    ElasticCableConstantGravity();

    libMesh::Point _gravity;
  };

} // end namespace GRINS

#endif /* GRINS_ELASTIC_CABLE_CONSTANT_GRAVITY_H */
