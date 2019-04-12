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


#ifndef GRINS_AVERAGE_NUSSELT_NUMBER_H
#define GRINS_AVERAGE_NUSSELT_NUMBER_H

// GRINS
#include "grins/qoi_base.h"
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  class PrimitiveTempFEVariables;

  class AverageNusseltNumber : public QoIBase
  {
  public:

    AverageNusseltNumber( const std::string& qoi_name );

    virtual ~AverageNusseltNumber();

    virtual QoIBase* clone() const;

    virtual bool assemble_on_interior() const;

    virtual bool assemble_on_sides() const;

    virtual void side_qoi( AssemblyContext& context,
                           const unsigned int qoi_index );

    virtual void side_qoi_derivative( AssemblyContext& context,
                                      const unsigned int qoi_index );

    virtual void init( const GetPot& input,
                       const MultiphysicsSystem& system,
                       unsigned int qoi_num );

    virtual void init_context( AssemblyContext& context );

  protected:

    void parse_thermal_conductivity( const GetPot& input );

    //! Thermal conductivity
    libMesh::Real _k;

    const PrimitiveTempFEVariables * _temp_vars;

    //! List of boundary ids for which we want to compute this QoI
    std::set<libMesh::boundary_id_type> _bc_ids;

    //! Scaling constant
    libMesh::Real _scaling;

  private:

    AverageNusseltNumber();

  };

  inline
  bool AverageNusseltNumber::assemble_on_interior() const
  {
    return false;
  }

  inline
  bool AverageNusseltNumber::assemble_on_sides() const
  {
    return true;
  }
}
#endif //GRINS_AVERAGE_NUSSELT_NUMBER_H
