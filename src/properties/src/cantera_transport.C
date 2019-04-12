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


#include "grins_config.h"

#ifdef GRINS_HAVE_CANTERA

// This class
#include "grins/cantera_transport.h"

// GRINS
#include "grins/cached_values.h"
#include "grins/cantera_mixture.h"

// libMesh
#include "libmesh/getpot.h"

// Cantera (with compiler warnings disabled)
#include "libmesh/ignore_warnings.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"
#include "libmesh/restore_warnings.h"

namespace GRINS
{

  CanteraTransport::CanteraTransport( CanteraMixture& mixture )
    : _cantera_gas( mixture.get_chemistry() ),
      _cantera_transport( mixture.get_transport() )
  {}

  libMesh::Real CanteraTransport::mu( const libMesh::Real& T,
                                      const libMesh::Real P,
                                      const std::vector<libMesh::Real>& Y )
  {
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    libMesh::Real mu = 0.0;

    {
      libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

      /*! \todo Need to make sure this will work in a threaded environment.
        Not sure if we will get thread lock here or not. */
      try
        {
          _cantera_gas.setState_TPY(T, P, &Y[0]);
          mu =  _cantera_transport.viscosity();
        }
      catch(Cantera::CanteraError)
        {
          Cantera::showErrors(std::cerr);
          libmesh_error();
        }

    }

    return mu;
  }

  libMesh::Real CanteraTransport::k( const libMesh::Real& T,
                                     const libMesh::Real P,
                                     const std::vector<libMesh::Real>& Y )
  {
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    libMesh::Real k = 0.0;

    {
      libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

      /*! \todo Need to make sure this will work in a threaded environment.
        Not sure if we will get thread lock here or not. */
      try
        {
          _cantera_gas.setState_TPY(T, P, &Y[0]);
          k =  _cantera_transport.thermalConductivity();
        }
      catch(Cantera::CanteraError)
        {
          Cantera::showErrors(std::cerr);
          libmesh_error();
        }

    }

    return k;
  }

  void CanteraTransport::mu_and_k_and_D( const libMesh::Real T,
                                         const libMesh::Real rho,
                                         const libMesh::Real /*cp*/,
                                         const std::vector<libMesh::Real>& Y,
                                         libMesh::Real& mu, libMesh::Real& k,
                                         std::vector<libMesh::Real>& D )
  {
    libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

    /*! \todo Need to make sure this will work in a threaded environment.
      Not sure if we will get thread lock here or not. */
    try
      {
        _cantera_gas.setState_TRY(T, rho, &Y[0]);

        mu =  _cantera_transport.viscosity();
        k =  _cantera_transport.thermalConductivity();
        _cantera_transport.getMixDiffCoeffsMass(&D[0]);
      }
    catch(Cantera::CanteraError)
      {
        Cantera::showErrors(std::cerr);
        libmesh_error();
      }
  }

} // namespace GRINS

#endif //GRINS_HAVE_CANTERA
