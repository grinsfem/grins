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
  template<typename KineticsThermoCurveFit>
  class AntiochMixture;

  //! Wrapper class for evaluating chemical kinetics using Antioch
  /*!
    This class is expected to be constructed *after* threads have been forked and will only
    live during the lifetime of the thread.
    By default, Antioch is working in SI units. Note that this documentation will always
    be built regardless if Antioch is included in the GRINS build or not. Check configure
    output to confirm that Antioch was included in the build.
  */
  template<typename KineticsThermoCurveFit>
  class AntiochKinetics
  {
  public:

    AntiochKinetics( const AntiochMixture<KineticsThermoCurveFit> & mixture );

    ~AntiochKinetics(){};

    void omega_dot( const libMesh::Real& T,
                    const libMesh::Real rho,
                    const std::vector<libMesh::Real>& mass_fractions,
                    std::vector<libMesh::Real>& omega_dot );

    void omega_dot( const Antioch::TempCache<libMesh::Real>& temp_cache,
                    const libMesh::Real rho,
                    const std::vector<libMesh::Real>& mass_fractions,
                    std::vector<libMesh::Real>& omega_dot );

  protected:

    const AntiochMixture<KineticsThermoCurveFit> & _antioch_mixture;

    Antioch::KineticsEvaluator<libMesh::Real> _antioch_kinetics;

    Antioch::NASAEvaluator<libMesh::Real,KineticsThermoCurveFit> _antioch_nasa_thermo;

  private:

    AntiochKinetics();

  };

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif //GRINS_ANTIOCH_KINETICS_H
