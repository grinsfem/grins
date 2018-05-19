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
#include "grins/average_nusselt_number.h"

// GRINS
#include "grins/multiphysics_sys.h"
#include "grins/assembly_context.h"
#include "grins/materials_parsing.h"
#include "grins/variable_warehouse.h"
#include "grins/variables_parsing.h"
#include "grins/single_variable.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"

namespace GRINS
{
  FlameSpeed::FlameSpeed( const std::string& qoi_name )
    : QoIBase(qoi_name)
  {
    return;
  }

  FlameSpeed::~FlameSpeed()
  {
    return;
  }

  QoIBase* FlameSpeed::clone() const
  {
    return new FlameSpeed( *this);
  }



  void FlameSpeed::init
  (const GetPot& input,
   const MultiphysicsSystem& /*system*/,
   unsigned int /*qoi_num*/ )
  {
    //grab the pressure value
    this->parse_Pressure(input);

    // Read boundary ids for which we want to compute, and make sure they are specified
    int num_bcs =  input.vector_variable_size("QoI/FlameSpeed/bc_ids");

    if( num_bcs <= 0 )
      {
        std::cerr << "Error: Must specify at least one boundary id to compute"
                  << " average Nusselt number." << std::endl
                  << "Found: " << num_bcs << std::endl;
        libmesh_error();
      }
    //insert the boundary id specified
    for( int i = 0; i < num_bcs; i++ )
      {
        _bc_ids.insert( input("QoI/FlameSpeed/bc_ids", -1, i ) );
      }

    //Set our Vaiables and number of species
    _temp_vars = &GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>(VariablesParsing::temp_variable_name(input,std::string("FlameSpeed"),VariablesParsing::QOI));
    _species_vars = &GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>(VariablesParsing::species_mass_frac_variable_name(input,std::string("FlameSpeed"),VariablesParsing::QOI));
    _mass_flux_vars = &GRINSPrivate::VariableWarehouse::get_variable_subclass<SingleVariable>(VariablesParsing::single_variable_name(input,std::string("FlameSpeed"),VariablesParsing::QOI));
    // #### might only need to use once, no need to set it _n_species = _species_vars.n_species();
  }






  ///Probaly NEED CHANGING OF SOME KIND
  void FlameSpeed::init_context( AssemblyContext& context )
  {
    libMesh::FEBase* T_fe;

    context.get_side_fe<libMesh::Real>(this->_temp_vars->T(), T_fe);

    T_fe->get_dphi();
    T_fe->get_JxW();

    return;
  }





  void FlameSpeed::side_qoi( AssemblyContext& context,
                                       const unsigned int qoi_index )
  {
    bool on_correct_side = false;

    for (std::set<libMesh::boundary_id_type>::const_iterator id =
           _bc_ids.begin(); id != _bc_ids.end(); id++ )
      if( context.has_side_boundary_id( (*id) ) )
        {
          on_correct_side = true;
          break;
        }

    if (!on_correct_side)
      return;

    libMesh::Number& qoi = context.get_qois()[qoi_index];
    libMesh::Real T = this->T(point,context);
    libMesh::Real p0 = this->_P0();
    libMesh::Real Mdot = this->M_dot(point,context);
    qoi += Mdot/(this->rho(T,p0, 500));
      } // quadrature loop
  }








  void FlameSpeed::parse_Pressure( const GetPot& input )
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
  }

} //namespace GRINS
