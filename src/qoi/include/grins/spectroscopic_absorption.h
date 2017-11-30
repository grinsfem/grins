//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_SPECTROSCOPIC_ABSORPTION_H
#define GRINS_SPECTROSCOPIC_ABSORPTION_H

// GRINS
#include "grins/integrated_function.h"
#include "grins/absorption_coeff.h"
#include "grins/single_variable.h"

namespace GRINS
{
  /*!
    QoI class for absorption calculations using the Beer-Lambert Law

    Relies upon the IntegratedFunction class for hooking into QoI infrastructure,
    and AbsorptionCoeff class for evaluating the <i>spectral absorption coefficient</i>, \f$ k_{\nu} \f$

    \f$ \frac{I_{\nu}}{I_{\nu}^0} = \exp\left\{- \int_0^L k_{\nu} dx\right\} \f$

    where \f$ \frac{I_{\nu}}{I_{\nu}^0} \f$ is the <i>spectral absorption</i>

    Expects all parameters given in standard SI units [m], [K], [Pa]
  */
  class SpectroscopicAbsorption : public IntegratedFunction<FEMFunctionAndDerivativeBase<libMesh::Real> >
  {
  public:

    /*!
      @param absorb AbsorptionCoeff object
    */
    SpectroscopicAbsorption(const GetPot & input,const std::string & qoi_name,SharedPtr<FEMFunctionAndDerivativeBase<libMesh::Real> > absorb);

    virtual QoIBase * clone() const;

    //! Override the QoIBase implementation to perform exp(-kv*L)
    virtual void parallel_op( const libMesh::Parallel::Communicator & communicator,
                              libMesh::Number & sys_qoi,
                              libMesh::Number & local_qoi );

  private:

    SpectroscopicAbsorption();

  };
}
#endif //GRINS_SPECTROSCOPIC_ABSORPTION_H
