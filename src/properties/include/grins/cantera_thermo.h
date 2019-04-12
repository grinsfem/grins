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


#ifndef GRINS_CANTERA_THERMO_H
#define GRINS_CANTERA_THERMO_H

#include "grins_config.h"

#ifdef GRINS_HAVE_CANTERA

// C++
#include <vector>

// libMesh
#include "libmesh/libmesh_common.h"

// libMesh forward declarations
class GetPot;

// Cantera forward declarations
namespace Cantera
{
  class IdealGasMix;
}

namespace GRINS
{
  // GRINS forward declarations
  class CanteraMixture;
  class CachedValues;

  //! Wrapper class for evaluating thermo properties using Cantera
  /*!
    This class is expected to be constructed *after* threads have been forked and will only
    live during the lifetime of the thread. Note that this documentation will always
    be built regardless if Cantera is included in the GRINS build or not. Check configure
    output to confirm that Cantera was included in the build if you wish to use it.
  */
  class CanteraThermodynamics
  {
  public:

    CanteraThermodynamics( CanteraMixture& mixture );
    ~CanteraThermodynamics(){};

    libMesh::Real cp( const libMesh::Real& T, const libMesh::Real P, const std::vector<libMesh::Real>& Y );
   
    void cp_s(const libMesh::Real& T,const libMesh::Real P, const std::vector<libMesh::Real>& Y , std::vector<libMesh::Real>& Cp_s);

    libMesh::Real cv( const libMesh::Real& T, const libMesh::Real P, const std::vector<libMesh::Real>& Y );

    libMesh::Real h( const libMesh::Real& T, unsigned int species );

  protected:

    CanteraMixture& _cantera_mixture;

    Cantera::IdealGasMix& _cantera_gas;

  private:

    CanteraThermodynamics();

  };

} // namespace GRINS

#endif //GRINS_HAVE_CANTERA

#endif //GRINS_CANTERA_THERMO_H
