//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_AVERAGE_NUSSELT_NUMBER_H
#define GRINS_AVERAGE_NUSSELT_NUMBER_H

// GRINS
#include "grins/qoi_base.h"
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  class AverageNusseltNumber : public QoIBase
  {
  public:

    AverageNusseltNumber( const GetPot& input );

    virtual ~AverageNusseltNumber();

    virtual libMesh::AutoPtr<libMesh::DifferentiableQoI> clone();

    virtual void side_qoi( libMesh::DiffContext& context, const libMesh::QoISet& qoi_indices );

    virtual void init( const GetPot& input, const MultiphysicsSystem& system );

  protected:

    virtual void read_input_options( const GetPot& input );

    //! Thermal conductivity
    libMesh::Real _k;

    //! Temperature variable index
    VariableIndex _T_var;

    //! List of boundary ids for which we want to compute this QoI
    std::set<libMesh::boundary_id_type> _bc_ids;

    //! Scaling constant
    libMesh::Real _scaling;

  private:

    AverageNusseltNumber();

  };
}
#endif //GRINS_AVERAGE_NUSSELT_NUMBER_H
