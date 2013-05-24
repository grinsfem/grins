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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_ANTIOCH_KINETICS_H
#define GRINS_ANTIOCH_KINETICS_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// C++
#include <vector>

// libMesh
#include "libmesh/libmesh_common.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/kinetics_evaluator.h"
#include "antioch/cea_evaluator.h"
namespace GRINS
{
  // GRINS forward declarations
  class AntiochMixture;

  class AntiochKinetics
  {
  public:

    AntiochKinetics( const AntiochMixture& mixture );
    ~AntiochKinetics();

    void omega_dot( const Antioch::TempCache<libMesh::Real>& temp_cache,
                    const libMesh::Real rho,
                    const libMesh::Real R_mix,
                    const std::vector<libMesh::Real>& mass_fractions,
                    std::vector<libMesh::Real>& omega_dot );

  protected:

    const AntiochMixture& _antioch_mixture;

    Antioch::KineticsEvaluator<libMesh::Real> _antioch_kinetics;

    Antioch::CEAEvaluator<libMesh::Real> _antioch_cea_thermo;

  private:

    AntiochKinetics();

  };

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif //GRINS_ANTIOCH_KINETICS_H
