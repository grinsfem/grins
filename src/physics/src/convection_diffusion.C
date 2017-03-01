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
#include "grins/convection_diffusion.h"

// GRINS
#include "grins/common.h"
#include "grins/assembly_context.h"
#include "grins/materials_parsing.h"
#include "grins/generic_ic_handler.h"
#include "grins/variable_warehouse.h"
#include "grins/single_variable.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  ConvectionDiffusion::ConvectionDiffusion( const GRINS::PhysicsName& physics_name,
                                            const GetPot& input )
    : Physics(physics_name,input),
      _v(3,libMesh::ParsedFunction<libMesh::Number>("0.0") ),
      _kappa("0.0"),
      _var(GRINSPrivate::VariableWarehouse::get_variable_subclass<SingleVariable>(VariablesParsing::single_variable_name(input,physics_name,VariablesParsing::PHYSICS)))
  {
    unsigned int n_v_comps = input.vector_variable_size("Physics/"+physics_name+"/velocity_field");

    for( unsigned int v = 0; v < n_v_comps; v++ )
      _v[v].reparse( input("Physics/"+physics_name+"/velocity_field", "0.0", v) );

    std::string material_name = MaterialsParsing::material_name(input,physics_name);

    this->set_parameter(this->_kappa, input,
                        "Materials/"+material_name+"/Diffusivity/value",
                        "DIE!");

    _ic_handler = new GenericICHandler(physics_name,input);

    this->check_var_subdomain_consistency(_var);
  }

  void ConvectionDiffusion::set_time_evolving_vars( libMesh::FEMSystem* system )
  {
    system->time_evolving(_var.var(),1);
  }

  void ConvectionDiffusion::init_context( AssemblyContext& context )
  {
    context.get_element_fe(_var.var())->get_JxW();
    context.get_element_fe(_var.var())->get_phi();
    context.get_element_fe(_var.var())->get_dphi();
    context.get_element_fe(_var.var())->get_xyz();
  }

  void ConvectionDiffusion::element_time_derivative
  ( bool compute_jacobian,
    AssemblyContext & context )
  {
    const unsigned int n_dofs = context.get_dof_indices(_var.var()).size();

    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_var.var())->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& phi =
      context.get_element_fe(_var.var())->get_phi();

    const std::vector<std::vector<libMesh::RealGradient> >& grad_phi =
      context.get_element_fe(_var.var())->get_dphi();

    const std::vector<libMesh::Point>& x =
      context.get_element_fe(_var.var())->get_xyz();

    libMesh::DenseSubVector<libMesh::Number> &F = context.get_elem_residual(_var.var());

    libMesh::DenseSubMatrix<libMesh::Number> &K = context.get_elem_jacobian(_var.var(), _var.var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Number c = context.interior_value(_var.var(), qp);

        libMesh::Gradient grad_c;
        context.interior_gradient(_var.var(), qp, grad_c);

        libMesh::RealGradient v_qp;
        v_qp(0) = this->_v[0](x[qp], context.get_time());
        v_qp(1) = this->_v[1](x[qp], context.get_time());
        v_qp(2) = this->_v[2](x[qp], context.get_time());

        libMesh::Real kappa_qp = this->_kappa(x[qp], context.time);

        for (unsigned int i=0; i != n_dofs; i++)
          {
            F(i) += JxW[qp]*( c*(v_qp*grad_phi[i][qp])
                              - kappa_qp*(grad_c*grad_phi[i][qp]) );

            if (compute_jacobian)
              {
                for (unsigned int j=0; j != n_dofs; j++)
                  {
                    K(i,j) += context.get_elem_solution_derivative()*
                      JxW[qp]*( phi[j][qp]*(v_qp*grad_phi[i][qp])
                                - kappa_qp*(grad_phi[j][qp]*grad_phi[i][qp]) );
                  }
              }
          }
      }
  }

  void ConvectionDiffusion::mass_residual
  ( bool compute_jacobian, AssemblyContext & context )
  {
    const unsigned int n_dofs = context.get_dof_indices(_var.var()).size();

    const std::vector<libMesh::Real> &JxW =
      context.get_element_fe(_var.var())->get_JxW();

    const std::vector<std::vector<libMesh::Real> >& phi =
      context.get_element_fe(_var.var())->get_phi();

    libMesh::DenseSubVector<libMesh::Number> &F = context.get_elem_residual(_var.var());

    libMesh::DenseSubMatrix<libMesh::Number> &M = context.get_elem_jacobian(_var.var(), _var.var());

    unsigned int n_qpoints = context.get_element_qrule().n_points();

    for (unsigned int qp=0; qp != n_qpoints; qp++)
      {
        libMesh::Real c_dot;
        context.interior_rate(_var.var(), qp, c_dot);

        for (unsigned int i=0; i != n_dofs; i++)
          {
            F(i) -= JxW[qp]*( c_dot*phi[i][qp] );

            if (compute_jacobian)
              {
                for (unsigned int j=0; j != n_dofs; j++)
                  {
                    M(i,j) -=
                      context.get_elem_solution_rate_derivative()*
                      JxW[qp]*( phi[j][qp]*phi[i][qp] );
                  }
              }
          }
      }
  }

} // end namespace GRINS
