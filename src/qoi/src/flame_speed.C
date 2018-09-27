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
#include "grins/flame_speed.h"


// GRINS
#include "grins/variable_warehouse.h"
#include "grins/chemistry_builder.h"
#include "grins/physics.h"

// libMesh
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"
#include "libmesh/elem.h"

#if GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

#if GRINS_HAVE_CANTERA
#include "grins/cantera_mixture.h"
#endif

namespace GRINS
{
  template< typename Chemistry>
  FlameSpeed<Chemistry>::FlameSpeed( const std::string& qoi_name)
    : QoIBase(qoi_name)
  {
    return;
  }

  template< typename Chemistry>
  FlameSpeed<Chemistry>::~FlameSpeed()
  {
    return;
  }

  template< typename Chemistry>
  QoIBase* FlameSpeed<Chemistry>::clone() const
  {
     return new FlameSpeed(*this);
  }


  template< typename Chemistry>
  void FlameSpeed<Chemistry>::init
  (const GetPot& input,
   const MultiphysicsSystem& /*system*/,
   unsigned int /*qoi_num*/ )
  {
    //grab the pressure value
    this->parse_Pressure(input);
    //grab the unburnt temperature value
    this->set_parameter
      ( _T_unburnt, input,"QoI/FlameSpeed/Unburnt_Temperature" , -1.0 );



    //figure out the chemistry
     std::string material = input("QoI/FlameSpeed/material","DIE!");
     _chemistry.reset( new Chemistry(input,material));

     //For a one dimensional premixed flame, the unburnt species values have to already be
     //specified for the physics, grabbing those values.
     Mass_Fractions.resize(this->_chemistry->n_species());
     for (unsigned int s = 0; s < this->_chemistry->n_species(); s++)
       {
         this->set_parameter(Mass_Fractions[s], input, "Physics/"+PhysicsNaming::od_premixed_flame()+"/"
                             +_chemistry->species_name(s), 0.0);
       }

    //Set our Vaiables
    _mass_flux_vars = &GRINSPrivate::VariableWarehouse::get_variable_subclass<SingleVariable>(VariablesParsing::single_variable_name(input,std::string("FlameSpeed"),VariablesParsing::QOI));

  }




  ///Probaly NEED CHANGING OF SOME KIND
  template< typename Chemistry>
  void FlameSpeed<Chemistry>::init_context( AssemblyContext& context )
  {
    libMesh::FEBase* M_fe;

    context.get_element_fe<libMesh::Real>(this->_mass_flux_vars->var(), M_fe);

    M_fe->get_dphi();
    M_fe->get_JxW();
    M_fe->get_xyz();

    return;
  }

  template< typename Chemistry>
  void FlameSpeed<Chemistry>::element_qoi( AssemblyContext & context,
                                           const unsigned qoi_index)
  {
    if(context.get_elem().contains_point(0.01) ) //if were at 1 cm off the boundary
      {

        libMesh::Number& qoi = context.get_qois()[qoi_index];
        libMesh::Real p0 = this->_P0;
        libMesh::Real Mdot = context.point_value(this->_mass_flux_vars->var(),0.01);


        libMesh::Real Rmix = _chemistry->R_mix(Mass_Fractions);

        qoi += Mdot/(this->rho(_T_unburnt,p0, Rmix));
      }
  }

  template< typename Chemistry>
  void FlameSpeed<Chemistry>::parse_Pressure( const GetPot& input)
  {
    //grab where the value is located in input file
    std::string material = input("QoI/FlameSpeed/material", "NoMaterial!");

    if( input.have_variable("Materials/"+material+"/ThermodynamicPressure/value") )
      {
        this->set_parameter
          ( _P0, input, "Materials/"+material+"/ThermodynamicPressure/value", -1.0 );
      }
    else
      libmesh_error_msg("ERROR: Could not find valid thermodynamicPressure value!");
  }//end parse pressure





  #if GRINS_HAVE_ANTIOCH
  template class FlameSpeed<AntiochChemistry>;
#endif

#if GRINS_HAVE_CANTERA
  template class FlameSpeed<CanteraMixture>;
  #endif
} //namespace GRINS
