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


// This class
#include "grins/spectroscopic_qoi_base.h"
#include "grins/absorption_coeff.h"
#include "grins/integrated_function.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"
#include "grins/materials_parsing.h"
#include "grins/single_variable.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_function_base.h"

namespace GRINS
{
  SpectroscopicQoIBase::SpectroscopicQoIBase( const std::shared_ptr<FEMFunctionAndDerivativeBase<libMesh::Real> > & absorb,
                                              const std::shared_ptr<RayfireMesh> & rayfire, const std::string & qoi_name, bool output_as_csv)
    : IntegratedFunction<FEMFunctionAndDerivativeBase<libMesh::Real> >(2 /* QGauss order */,absorb,rayfire,qoi_name),
      _output_as_csv(output_as_csv)
  {}

  void SpectroscopicQoIBase::output_qoi(std::ostream & out) const
  {
    if (_output_as_csv)
      {
        const AbsorptionCoeffBase & abs = libMesh::cast_ref<const AbsorptionCoeffBase &>(this->get_function());
        libMesh::Real nu = abs.get_wavenumber();

        out << std::setprecision(16)
            << std::scientific
            << nu << ","
            << _qoi_value << std::endl;
      }
    else
      {
        QoIBase::output_qoi(out);
      }
  }

} //namespace GRINS
