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
  template< typename Chemistry>
  FlameSpeed<Chemistry>::FlameSpeed( const std::string& qoi_name, std::unique_ptr<Chemistry> & chem)
    : QoIBase(qoi_name),
      _chemistry(chem.release())
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
    libmesh_not_implemented() ;
  }


  template< typename Chemistry>
  void FlameSpeed<Chemistry>::init
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
    _mass_flux_vars = &GRINSPrivate::VariableWarehouse::get_variable_subclass<SingleVariable>(VariablesParsing::single_variable_name(input,std::string("FlameSpeed"),VariablesParsing::QOI));
    _species_vars= &GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>(VariablesParsing::species_mass_frac_variable_name(input,std::string("FlameSpeed"),VariablesParsing::QOI));
    // #### might only need to use once, no need to set it _n_species = _species_vars.n_species();
  }






  ///Probaly NEED CHANGING OF SOME KIND
  template< typename Chemistry>
  void FlameSpeed<Chemistry>::init_context( AssemblyContext& context )
  {
    libMesh::FEBase* T_fe;

    context.get_side_fe<libMesh::Real>(this->_temp_vars->T(), T_fe);

    T_fe->get_dphi();
    T_fe->get_JxW();

    return;
  }


  //###############################################################################################################

  template< typename Chemistry>
  void FlameSpeed<Chemistry>::side_qoi( AssemblyContext& context,
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
    libMesh::Real T = context.side_value(this->_temp_vars->T(),0);
    libMesh::Real p0 = this->_P0;
    libMesh::Real Mdot = context.side_value(this->_mass_flux_vars->var(),0);

    std::vector<libMesh::Real> mass_fractions;
    mass_fractions.resize(this->_species_vars->n_species());
    for (unsigned int s=0; s < this->_species_vars->n_species(); s++)
      {
	mass_fractions[s] = context.side_value(this->_species_vars->species(s),0);
      }
    libMesh::Real Rmix = _chemistry->R_mix(mass_fractions);

    qoi += Mdot/(this->rho(T,p0, Rmix));
  } //end side_qoi


  //###########################################################################################################
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

} //namespace GRINS
