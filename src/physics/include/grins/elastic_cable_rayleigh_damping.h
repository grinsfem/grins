//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_ELASTIC_CABLE_RAYLEIGH_DAMPING_H
#define GRINS_ELASTIC_CABLE_RAYLEIGH_DAMPING_H

#include "grins/elastic_cable_base.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  class ElasticCableRayleighDamping : public ElasticCableBase<StressStrainLaw>
  {
  public:

    ElasticCableRayleighDamping( const PhysicsName& physics_name,
                                 const GetPot& input,
                                 bool is_compressible );

    virtual ~ElasticCableRayleighDamping(){};

    //! Error out if using libMesh::FirstOrderUnsteadySolver
    virtual void auxiliary_init( MultiphysicsSystem & system );

    //! Time dependent part(s) of physics for element interiors
    virtual void damping_residual( bool compute_jacobian,
                                   AssemblyContext & context );

  protected:

    libMesh::Real _lambda_factor;
    libMesh::Real _mu_factor;

  private:

    ElasticCableRayleighDamping();

  };

} // end namespace GRINS

#endif // GRINS_ELASTIC_CABLE_RAYLEIGH_DAMPING_H
