//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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


#ifndef ADIABIATIC_FLAME_TEMPERATURE_H
#define ADIABIATIC_FLAME_TEMPERATURE_H

// GRINS
#include "grins/qoi_base.h"

#include "grins/multiphysics_sys.h"
#include "grins/single_variable.h"
#include "grins/multicomponent_variable.h"

namespace GRINS
{
  class PrimitiveTempFEVariables;
  class AdiabiaticFlameTemperature : public QoIBase
  {
  public:

    AdiabiaticFlameTemperature(const std::string& qoi_name);

    virtual ~AdiabiaticFlameTemperature();

    virtual QoIBase* clone() const;

    virtual bool assemble_on_interior() const;

    virtual bool assemble_on_sides() const;

    virtual void element_qoi( AssemblyContext& context,
                              const unsigned int qoi_index);



    virtual void init( const GetPot& input,
                       const MultiphysicsSystem& system,
                       unsigned int qoi_num );

    virtual void init_context( AssemblyContext& context );

    libMesh::Real T( const libMesh::Point& p, const AssemblyContext& c ) const;
  protected:
    const PrimitiveTempFEVariables * _temp_vars;

  private:

    AdiabiaticFlameTemperature();

  };

  inline
    libMesh::Real AdiabiaticFlameTemperature::T( const libMesh::Point& p,
                                                 const AssemblyContext& c ) const
  { return c.point_value(_temp_vars->var(),p); }



    inline
      bool AdiabiaticFlameTemperature::assemble_on_interior() const
    {
      return true;
    }

    inline
      bool AdiabiaticFlameTemperature::assemble_on_sides() const
    {
      return false;
    }


}
#endif //GRINS_ADIABIATIC_FLAME_TEMPERATURE_H
