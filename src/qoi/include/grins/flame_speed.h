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
#include "grins_config.h"
#include "grins/grins_enums.h"

#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"
#include "grins/materials_parsing.h"
#include "grins/variable_warehouse.h"
#include "grins/variables_parsing.h"
#include "grins/single_variable.h"
#include "grins/multi_component_vector_variable.h"
#include "grins/multicomponent_variable.h"

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  class PrimitiveTempFEVariables;
  template<typename Chemistry>
    class FlameSpeed : public QoIBase
  {
  public:

    FlameSpeed(const std::string& qoi_name, std::unique_ptr<Chemistry> & chem);
    
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
    void mass_fractions( const libMesh::Point& p, const AssemblyContext& c,
                         std::vector<libMesh::Real>& mass_fracs ) const;



  protected:

    void parse_Pressure( const GetPot& input);

    const PrimitiveTempFEVariables * _temp_vars;
    const SingleVariable * _mass_flux_vars;
    const SpeciesMassFractionsVariable * _species_vars;
    
    libMesh::Real rho( libMesh::Real T, libMesh::Real p0, libMesh::Real R_mix) const;

    //! List of boundary ids for which we want to compute this QoI
    std::set<libMesh::boundary_id_type> _bc_ids;

 
    std::unique_ptr<Chemistry> _chemistry;

    libMesh::Real _P0;

  private:

    FlameSpeed();

  };

  template< typename Chemistry>
    inline
    libMesh::Real FlameSpeed<Chemistry>::rho( libMesh::Real T,
				   libMesh::Real p0,
				   libMesh::Real R_mix) const
  {
    libMesh::Real value = 0;
    value = p0/(R_mix*T);
    return value;
  }
  
  template< typename Chemistry>
    inline
    libMesh::Real FlameSpeed<Chemistry>::M_dot( const libMesh::Point& p,
				     const AssemblyContext& c ) const
  { return c.point_value(_mass_flux_vars->var(),p); }
  
  template< typename Chemistry>
    inline
    libMesh::Real FlameSpeed<Chemistry>::T( const libMesh::Point& p,
				 const AssemblyContext& c ) const
  { return c.point_value(_temp_vars->T(),p); }
  template< typename Chemistry>
    inline
    bool FlameSpeed<Chemistry>::assemble_on_interior() const
  {
    return false;
  }
  template< typename Chemistry>
  inline
    bool FlameSpeed<Chemistry>::assemble_on_sides() const
  {
   return true;
  }
  
    
}
#endif //GRINS_Flame_Speed_H
