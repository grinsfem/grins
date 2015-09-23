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


#ifndef GRINS_WEIGHTED_FLUX_QOI_H
#define GRINS_WEIGHTED_FLUX_QOI_H

// GRINS
#include "grins/qoi_base.h"
#include "grins/variable_name_defaults.h"

// libMesh
#include "libmesh/auto_ptr.h"
#include "libmesh/function_base.h"

// C++
#include <set>
#include <string>
#include <vector>

namespace GRINS
{
  class WeightedFluxQoI : public QoIBase
  {
  public:

    WeightedFluxQoI( const std::string& qoi_name );

    virtual ~WeightedFluxQoI();

    virtual QoIBase* clone() const;

    virtual bool assemble_on_interior() const;

    virtual bool assemble_on_sides() const;

    virtual void init( const GetPot& input,
                       const MultiphysicsSystem& system,
                       unsigned int qoi_num );

  protected:

    //! List of variables which we want to contribute to fluxes
    std::vector<std::string> _var_names;

    //! List of boundary ids on which we want to compute this QoI
    std::set<libMesh::boundary_id_type> _bc_ids;

    //! Boundary condition functor
    libMesh::AutoPtr<libMesh::FunctionBase<libMesh::Number> >
      _adjoint_bc;

    //! Manual copy constructor due to the AutoPtr
    WeightedFluxQoI(const WeightedFluxQoI& original);

  private:

    WeightedFluxQoI();

  };

  inline
  bool WeightedFluxQoI::assemble_on_interior() const
  {
    // Although we are technically evaluating a lift-function-weighted
    // residual on the interior, our evaluation gets done by the
    // FEMSystem
    return false;
  }

  inline
  bool WeightedFluxQoI::assemble_on_sides() const
  {
    // Although our QoI is a boundary integral, it doesn't get
    // evaluated directly as such; FEMSystem uses a superconvergent
    // flux calculation based on the adjoint BC we added.
    return false;
  }
}
#endif //GRINS_WEIGHTED_FLUX_QOI_H
