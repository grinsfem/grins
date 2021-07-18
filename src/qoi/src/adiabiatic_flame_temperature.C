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


// This class
#include "grins/adiabiatic_flame_temperature.h"


// GRINS
#include "grins/variable_warehouse.h"
#include "grins/chemistry_builder.h"
#include "grins/physics.h"

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"
#include "libmesh/elem.h"


namespace GRINS
{
  AdiabiaticFlameTemperature::AdiabiaticFlameTemperature( const std::string& qoi_name)
    : QoIBase(qoi_name)
  {
    return;
  }

  AdiabiaticFlameTemperature::~AdiabiaticFlameTemperature()
  {
    return;
  }

  QoIBase* AdiabiaticFlameTemperature::clone() const
  {
     return new AdiabiaticFlameTemperature(*this);
  }


  void AdiabiaticFlameTemperature::init
  (const GetPot& input,
   const MultiphysicsSystem& /*system*/,
   unsigned int /*qoi_num*/ )
  {
    _temp_vars = &GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>(VariablesParsing::temp_variable_name(input,std::string("AdiabiaticFlameTemperature"),VariablesParsing::QOI));
    this->set_parameter
      ( _T_point, input, "QoI/AdiabiaticFlameTemperature/point", .01 );
  }




  ///Probaly NEED CHANGING OF SOME KIND
  void AdiabiaticFlameTemperature::init_context( AssemblyContext& context )
  {
    libMesh::FEBase* T_fe;

    context.get_element_fe<libMesh::Real>(this->_temp_vars->var(), T_fe);

    T_fe->get_dphi();
    T_fe->get_JxW();
    T_fe->get_xyz();

    return;
  }

  void AdiabiaticFlameTemperature::element_qoi( AssemblyContext & context,
                                           const unsigned qoi_index)
  {
    if(context.get_elem().contains_point(this->_T_point) ) //if were at 1 cm off the boundary
      {

        libMesh::Number& qoi = context.get_qois()[qoi_index];
        libMesh::Real T = context.point_value(this->_temp_vars->var(),0.01);
        qoi += T;
      }
  }


} //namespace GRINS
