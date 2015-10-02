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


#ifndef GRINS_BOUSSINESQ_BUOYANCY_BASE_H
#define GRINS_BOUSSINESQ_BUOYANCY_BASE_H

// GRINS
#include "grins/physics.h"
#include "grins/primitive_flow_fe_variables.h"
#include "grins/primitive_temp_fe_variables.h"

// libMesh
#include "libmesh/point.h"

namespace GRINS
{  
  class BoussinesqBuoyancyBase : public Physics
  {
  public:
    
    BoussinesqBuoyancyBase( const std::string& physics_name, const GetPot& input );

    ~BoussinesqBuoyancyBase();

    //! Initialization of BoussinesqBuoyancy variables
    virtual void init_variables( libMesh::FEMSystem* system );

  protected:

    //! Helper function for parsing/maintaing backward compatibility
    void read_property( const GetPot& input,
                        const std::string& old_option,
                        const std::string& property,
                        libMesh::Real& value );

    //! Helper function for parsing/maintaing backward compatibility
    void duplicate_input_test( const GetPot& input,
                               const std::string& option1,
                               const std::string& option2 );

    //! Helper function for parsing/maintaing backward compatibility
    void no_input_warning( const GetPot& input,
                           const std::string& old_option,
                           const std::string& material,
                           const std::string& property );

    PrimitiveFlowFEVariables _flow_vars;
    PrimitiveTempFEVariables _temp_vars;

    //! \f$ \rho = \f$ density
    libMesh::Number _rho;

    //! \f$ T_0 = \f$ reference temperature 
    libMesh::Number _T_ref;

    //! \f$ \beta_T = \f$ coefficient of thermal expansion
    libMesh::Number _beta_T;

    //! Gravitational vector
    libMesh::Point _g;

     //! Physical dimension of problem
    unsigned int _dim;

  private:

    BoussinesqBuoyancyBase();

  };

} // end namespace GRINS
#endif // GRINS_BOUSSINESQ_BUOYANCY_BASE_H
