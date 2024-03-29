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


#ifndef GRINS_SPECTROSCOPIC_TRANSMISSION_H
#define GRINS_SPECTROSCOPIC_TRANSMISSION_H

// GRINS
#include "grins/spectroscopic_qoi_base.h"
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

    where \f$ \frac{I_{\nu}}{I_{\nu}^0} \f$ is the <i>spectral transmission</i>

    Expects all parameters given in standard SI units [m], [K], [Pa]
  */
  class SpectroscopicTransmission : public SpectroscopicQoIBase
  {
  public:

    /*!
      @param absorb AbsorptionCoeff object
      @param output_as_csv Flag for whether we should output QoI value in wavenumber,absorption CSV format
        or in the normal QoIBase way
    */
    SpectroscopicTransmission(const std::shared_ptr<FEMFunctionAndDerivativeBase<libMesh::Real> > & absorb,
                               const std::shared_ptr<RayfireMesh> & rayfire, const std::string & qoi_name, bool output_as_csv);

    virtual QoIBase * clone() const override;

    //! Override the QoIBase implementation to perform exp(-kv*L)
    virtual void parallel_op( const libMesh::Parallel::Communicator & communicator,
                              libMesh::Number & sys_qoi,
                              libMesh::Number & local_qoi ) override;

    //! Override DifferentiableQoI's empty implementation to add chain rule (QoI is exponential)
    virtual void finalize_derivative(libMesh::NumericVector<libMesh::Number> & derivatives, std::size_t qoi_index) override;



  };
}
#endif //GRINS_SPECTROSCOPIC_TRANSMISSION_H
