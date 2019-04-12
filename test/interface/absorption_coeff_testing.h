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


#ifndef GRINS_ABSORPTION_COEFF_TESTING_H
#define GRINS_ABSORPTION_COEFF_TESTING_H

// GRINS
#include "grins/absorption_coeff.h"

namespace GRINSTesting
{
  template<typename Chemistry>
  class AbsorptionCoeffTesting : public GRINS::AbsorptionCoeff<Chemistry>
  {
  public:
    /*!
      This class is intended solely for use in unit testing the AbsorptionCoeff class.
    */
    AbsorptionCoeffTesting( std::shared_ptr<Chemistry> & chem, std::shared_ptr<GRINS::HITRAN> & hitran,
                            libMesh::Real nu_min, libMesh::Real nu_max,
                            libMesh::Real desired_nu, const std::string & species,
                            libMesh::Real thermo_pressure);

    friend class SpectroscopicTestBase;
  };

  template<typename Chemistry>
  AbsorptionCoeffTesting<Chemistry>::AbsorptionCoeffTesting(std::shared_ptr<Chemistry> & chem,
                                                            std::shared_ptr<GRINS::HITRAN> & hitran,
                                                            libMesh::Real nu_min, libMesh::Real nu_max,
                                                            libMesh::Real desired_nu, const std::string & species,
                                                            libMesh::Real thermo_pressure)
    : GRINS::AbsorptionCoeff<Chemistry>(chem,hitran,nu_min,nu_max,desired_nu,species,thermo_pressure)
  {}
}

#endif //GRINS_ABSORPTION_COEFF_TESTING_H
