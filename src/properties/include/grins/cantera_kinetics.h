//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_CANTERA_KINETICS_H
#define GRINS_CANTERA_KINETICS_H

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
  class CachedValues;
  class CanteraMixture;

  class CanteraKinetics
  {
  public:

    CanteraKinetics( CanteraMixture& mixture );
    ~CanteraKinetics();

    void omega_dot( const CachedValues& cache, unsigned int qp,
		    std::vector<libMesh::Real>& omega_dot ) const;

    void omega_dot_TPY( const libMesh::Real T, const libMesh::Real P,
                        const std::vector<libMesh::Real>& mass_fractions,
                        std::vector<libMesh::Real>& omega_dot ) const;

    void omega_dot_TRY( const libMesh::Real& T, const libMesh::Real rho,
                        const std::vector<libMesh::Real>& mass_fractions,
                        std::vector<libMesh::Real>& omega_dot ) const;

  protected:

    Cantera::IdealGasMix& _cantera_gas;

  private:

    CanteraKinetics();

  };

} // namespace GRINS

#endif // GRINS_HAVE_CANTERA

#endif //GRINS_CANTERA_KINETICS_H
