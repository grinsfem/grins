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
  AverageNusseltNumber::AverageNusseltNumber( const std::string& qoi_name )
    : QoIBase(qoi_name)
  {
    return;
  }

  AverageNusseltNumber::~AverageNusseltNumber()
  {
    return;
  }

  QoIBase* AverageNusseltNumber::clone() const
  {
    AverageNusseltNumber *returnval = new AverageNusseltNumber( *this );
    returnval->move_parameter(_k, returnval->_k);
    returnval->move_parameter(_scaling, returnval->_scaling);
    return returnval;
  }

  void AverageNusseltNumber::init
  (const GetPot& input,
   const MultiphysicsSystem& /*system*/,
   unsigned int /*qoi_num*/ )
  {
    this->parse_thermal_conductivity(input);

    this->set_parameter
      ( _scaling, input, "QoI/NusseltNumber/scaling", 1.0 );

    if( this->_k < 0.0 )
      {
        std::cerr << "Error: thermal conductivity for AverageNusseltNumber must be positive." << std::endl
                  << "Found k = " << _k << std::endl;
        libmesh_error();
      }

    // Read boundary ids for which we want to compute
    int num_bcs =  input.vector_variable_size("QoI/NusseltNumber/bc_ids");

    if( num_bcs <= 0 )
      {
        std::cerr << "Error: Must specify at least one boundary id to compute"
                  << " average Nusselt number." << std::endl
                  << "Found: " << num_bcs << std::endl;
        libmesh_error();
      }

    for( int i = 0; i < num_bcs; i++ )
      {
        _bc_ids.insert( input("QoI/NusseltNumber/bc_ids", -1, i ) );
      }

    _temp_vars = &GRINSPrivate::VariableWarehouse::get_variable_subclass<PrimitiveTempFEVariables>(VariablesParsing::temp_variable_name(input,std::string("NusseltNumber"),VariablesParsing::QOI));
  }

  void AverageNusseltNumber::init_context( AssemblyContext& context )
  {
    libMesh::FEBase* T_fe;

    context.get_side_fe<libMesh::Real>(this->_temp_vars->T(), T_fe);

    T_fe->get_dphi();
    T_fe->get_JxW();

    return;
  }

  void AverageNusseltNumber::side_qoi( AssemblyContext& context,
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

    libMesh::FEBase* side_fe;
    context.get_side_fe<libMesh::Real>(this->_temp_vars->T(), side_fe);

    const std::vector<libMesh::Real> &JxW = side_fe->get_JxW();

    const std::vector<libMesh::Point>& normals = side_fe->get_normals();

    unsigned int n_qpoints = context.get_side_qrule().n_points();

    libMesh::Number& qoi = context.get_qois()[qoi_index];

    // Loop over quadrature points

    for (unsigned int qp = 0; qp != n_qpoints; qp++)
      {
        // Get the solution value at the quadrature point
        libMesh::Gradient grad_T = 0.0;
        context.side_gradient(this->_temp_vars->T(), qp, grad_T);

        // Update the elemental increment dR for each qp
        qoi += (this->_scaling)*(this->_k)*(grad_T*normals[qp])*JxW[qp];

      } // quadrature loop
  }

  void AverageNusseltNumber::side_qoi_derivative( AssemblyContext& context,
                                                  const unsigned int qoi_index )
  {

    for( std::set<libMesh::boundary_id_type>::const_iterator id = _bc_ids.begin();
         id != _bc_ids.end(); id++ )
      {
        if( context.has_side_boundary_id( (*id) ) )
          {
            libMesh::FEBase* T_side_fe;
            context.get_side_fe<libMesh::Real>(this->_temp_vars->T(), T_side_fe);

            const std::vector<libMesh::Real> &JxW = T_side_fe->get_JxW();

            const std::vector<libMesh::Point>& normals = T_side_fe->get_normals();

            unsigned int n_qpoints = context.get_side_qrule().n_points();

            const unsigned int n_T_dofs = context.get_dof_indices(this->_temp_vars->T()).size();

            const std::vector<std::vector<libMesh::Gradient> >& T_gradphi = T_side_fe->get_dphi();

            libMesh::DenseSubVector<libMesh::Number>& dQ_dT =
              context.get_qoi_derivatives(qoi_index, this->_temp_vars->T());

            // Loop over quadrature points
            for (unsigned int qp = 0; qp != n_qpoints; qp++)
              {
                // Get the solution value at the quadrature point
                libMesh::Gradient grad_T = 0.0;
                context.side_gradient(this->_temp_vars->T(), qp, grad_T);

                // Update the elemental increment dR for each qp
                //qoi += (this->_scaling)*(this->_k)*(grad_T*normals[qp])*JxW[qp];

                for( unsigned int i = 0; i != n_T_dofs; i++ )
                  {
                    dQ_dT(i) += _scaling*_k*T_gradphi[i][qp]*normals[qp]*JxW[qp];
                  }

              } // quadrature loop

          } // end check on boundary id

      }

    return;
  }

  void AverageNusseltNumber::parse_thermal_conductivity( const GetPot& input )
  {
    std::string material = input("QoI/NusseltNumber/material", "NoMaterial!");

    MaterialsParsing::duplicate_input_test(input,
                                           "Materials/"+material+"/ThermalConductivity/value",
                                           "QoI/NusseltNumber/thermal_conductivity");

    // Parse the old version
    if( input.have_variable("QoI/NusseltNumber/thermal_conductivity") )
      {
        MaterialsParsing::dep_input_warning( "QoI/NusseltNumber/thermal_conductivity",
                                             "ThermalConductivity/value" );

        this->set_parameter
          ( _k, input, "QoI/NusseltNumber/thermal_conductivity", -1.0 );
      }
    // Parse new version
    else if( input.have_variable("Materials/"+material+"/ThermalConductivity/value") )
      {
        // This is currently only valid for constant thermal conductivity models
        if( input("Materials/"+material+"/ThermalConductivity/model", "DIE!")
            != std::string("constant") )
          {
            libmesh_error_msg("ERROR: Only constant ThermalConductivity model supported in NusseltNumber!");
          }

        this->set_parameter
          ( _k, input, "Materials/"+material+"/ThermalConductivity/value", -1.0 );
      }
    else
      libmesh_error_msg("ERROR: Could not find valid thermal conducitivity value!");

    // thermal conducivity should be positive
    if( _k <= 0.0 )
      libmesh_error_msg("ERROR: Detected non-positive thermal conductivity!");
  }

} //namespace GRINS
