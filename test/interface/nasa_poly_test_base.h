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

#ifndef GRINS_NASA_POLY_TEST_BASE_H
#define GRINS_NASA_POLY_TEST_BASE_H

namespace GRINSTesting
{
  class NASA7TestBase
  {
  public:

    libMesh::Real cp_R_exact( libMesh::Real T,
                              libMesh::Real a0, libMesh::Real a1, libMesh::Real a2,
                              libMesh::Real a3, libMesh::Real a4 )
    {
      return a0 + a1*T + a2*T*T + a3*T*T*T + a4*(T*T*T*T);
    }

    libMesh::Real h_RT_exact( libMesh::Real T,
                              libMesh::Real a0, libMesh::Real a1, libMesh::Real a2,
                              libMesh::Real a3, libMesh::Real a4, libMesh::Real a5 )
    {
      return a0 + a1/2.0*T + a2/3.0*T*T + a3/4.0*T*T*T + a4/5.0*(T*T*T*T) + a5/T;
    }

    libMesh::Real s_R_exact( libMesh::Real T,
                             libMesh::Real a0, libMesh::Real a1, libMesh::Real a2,
                             libMesh::Real a3, libMesh::Real a4, libMesh::Real a6 )
    {
      return a0*std::log(T) + a1*T + a2/2.0*T*T + a3/3.0*T*T*T + a4/4.0*(T*T*T*T) + a6;
    }
  };

  class NASA9TestBase
  {
  public:

    libMesh::Real cp_R_exact( libMesh::Real T,
                              libMesh::Real a0, libMesh::Real a1, libMesh::Real a2,
                              libMesh::Real a3, libMesh::Real a4, libMesh::Real a5, libMesh::Real a6 )
    {
      return a0/(T*T) + a1/T + a2 + a3*T + a4*(T*T) + a5*(T*T*T) + a6*(T*T*T*T);
    }

    libMesh::Real h_RT_exact( libMesh::Real T,
                              libMesh::Real a0, libMesh::Real a1, libMesh::Real a2, libMesh::Real a3,
                              libMesh::Real a4, libMesh::Real a5, libMesh::Real a6, libMesh::Real a7 )
    {
      return -a0/(T*T) + a1*std::log(T)/T + a2 + a3*T/2.0 + a4*(T*T)/3.0
        + a5*(T*T*T)/4.0 + a6*(T*T*T*T)/5.0 + a7/T;
    }

    libMesh::Real s_R_exact( libMesh::Real T,
                             libMesh::Real a0, libMesh::Real a1, libMesh::Real a2, libMesh::Real a3,
                             libMesh::Real a4, libMesh::Real a5, libMesh::Real a6, libMesh::Real a8 )
    {
      return -a0/(2.*T*T) - a1/T + a2*std::log(T) + a3*T + a4*(T*T)/2.0
        + a5*(T*T*T)/3.0 + a6*(T*T*T*T)/4.0 + a8;
    }
  };

} // end namespace GRINSTesting

#endif // GRINS_NASA_POLY_TEST_BASE_H
