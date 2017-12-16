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


// This class
#include "grins/stab_helper.h"

// GRINS
#include "grins/assembly_context.h"

namespace GRINS
{

  StabilizationHelper::StabilizationHelper
  ( const std::string & helper_name ) :
    ParameterUser(helper_name)
  {
    return;
  }

  StabilizationHelper::~StabilizationHelper()
  {
    return;
  }

  libMesh::RealGradient StabilizationHelper::compute_g( libMesh::FEBase* fe,
                                                        AssemblyContext& c,
                                                        unsigned int qp ) const
  {
    libMesh::RealGradient g( fe->get_dxidx()[qp] );

    if( c.get_dim() > 1 )
      {
        g(0) += fe->get_detadx()[qp];
        g(1) = fe->get_dxidy()[qp] + fe->get_detady()[qp];
      }

    if( c.get_dim() == 3 )
      {
        g(0) += fe->get_dzetadx()[qp];
        g(1) += fe->get_dzetady()[qp];
        g(2) = fe->get_dxidz()[qp] + fe->get_detadz()[qp] + fe->get_dzetadz()[qp];
      }

    return g;
  }

  libMesh::RealTensor StabilizationHelper::compute_G( libMesh::FEBase* fe,
                                                      AssemblyContext& c,
                                                      unsigned int qp ) const
  {
    libMesh::Real dxidx = fe->get_dxidx()[qp];

    libMesh::RealTensor G( dxidx*dxidx );

    if( c.get_dim() > 1 )
      {
        libMesh::Real dxidy = fe->get_dxidy()[qp];

        libMesh::Real detadx = fe->get_detadx()[qp];
        libMesh::Real detady = fe->get_detady()[qp];

        G(0,0) += detadx*detadx;
        G(0,1) = G(1,0) = dxidx*dxidy + detadx*detady;
        G(1,1) = dxidy*dxidy + detady*detady;

        if( c.get_dim() == 3 )
          {
            libMesh::Real dxidz = fe->get_dxidz()[qp];

            libMesh::Real detadz = fe->get_detadz()[qp];

            libMesh::Real dzetadx = fe->get_dzetadx()[qp];
            libMesh::Real dzetady = fe->get_dzetady()[qp];
            libMesh::Real dzetadz = fe->get_dzetadz()[qp];

            G(0,0) += dzetadx*dzetadx;
            G(0,1) += dzetadx*dzetady;
            G(0,2) = dxidx*dxidz + detadx*detadz + dzetadx*dzetadz;
            G(1,0) += dzetady*dzetadx;
            G(1,1) += dzetady*dzetady;
            G(1,2) = dxidy*dxidz + detady*detadz + dzetady*dzetadz;
            G(2,0) = dxidz*dxidx + detadz*detadx + dzetadz*dzetadx;
            G(2,1) = dxidz*dxidy + detadz*detady + dzetadz*dzetady;
            G(2,2) = dxidz*dxidz + detadz*detadz + dzetadz*dzetadz;
          }
      }

    return G;
  }

} // namespace GRINS
