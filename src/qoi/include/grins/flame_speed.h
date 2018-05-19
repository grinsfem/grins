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


#ifndef GRINS_FLAME_SPEED_H
#define GRINS_FLAME_SPEED_H

// GRINS
#include "grins/qoi_base.h"
#include "grins/variable_name_defaults.h"

namespace GRINS
{
  class PrimitiveTempFEVariables;

  class FlameSpeed : public QoiBase
  {
  public:

    FlameSpeed(const std::string& qoi_name);

    virtual ~FlameSpeed();

    virtual QoIBase* clone() const;

    virtual bool assemble_on_interior() const;

    virtual bool assemble_on_sides() const;

    virtual void side_qoi( AssemblyContext& context,
                           const unsigned int qoi_index );

    virtual void init( const GetPot& input,
                       const MultiphysicsSystem& system,
                       unsigned int qoi_num );

    virtual void init_context( AssemblyContext& context );

    libMesh::Real T( const libMesh::Point& p, const AssemblyContext& c ) const;
    libMesh::Real M_dot( const libMesh::Point& p, const AssemblyContext& c ) const;

    const Mixture & gas_mixture() const;

  protected:

    const PrimitiveTempFEVariables * _temp_vars;
    const SpeciesMassFractionsVariable * _species_vars;
    const SingleVariable * _mass_flux_vars;

    void parse_Pressure( const GetPot& input)

    const PrimitiveTempFEVariables * _temp_vars;

    //! List of boundary ids for which we want to compute this QoI
    libMesh::boundary_id_type _bc_ids;

    //! Scaling constant
    libMesh::Real _P0

  private:

    FlameSpeed();

  };

 

  inline
    libMesh::Real FlameSpeed::rho( libMesh::Real T,
				   libMesh::Real p0,
				   libMesh::Real R_mix) const
  {
    libMesh::Real value = 0;
    value = p0/(R_mix*T);
    return value;
  }
  
  inline
    libMesh::Real FlameSpeed::M_dot( const libMesh::Point& p,
							     const AssemblyContext& c ) const
    { return c.point_value(_mass_flux_vars.var(),p); }

  inline
    libMesh::Real FlameSpeed::T( const libMesh::Point& p,
							 const AssemblyContext& c ) const
    { return c.point_value(_temp_vars.T(),p); }
  
  inline
    bool FlameSpeed::assemble_on_interior() const
    {
      return false;
    }
  
  inline
    bool FlameSpeed::assemble_on_sides() const
    {
      return true;
    }
  
}
#endif //GRINS_Flame_Speed_H
