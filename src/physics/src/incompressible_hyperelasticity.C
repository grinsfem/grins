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
#include "grins/incompressible_hyperelasticity.h"

// GRINS
#include "grins/variables_parsing.h"
#include "grins/variable_warehouse.h"
#include "grins/multiphysics_sys.h"
#include "grins/cartesian_hyperelasticity.h"
#include "grins/incompressible_hyperelasticity_weak_form.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{
  template<unsigned int Dim,typename StrainEnergy>
  IncompressibleHyperelasticity<Dim,StrainEnergy>::IncompressibleHyperelasticity( const PhysicsName & physics_name,
                                                                                  const GetPot & input )
    :  HyperelasticityBase<Dim,StrainEnergy>(physics_name,physics_name,input),
    _press_var(GRINSPrivate::VariableWarehouse::get_variable_subclass<PressureFEVariable>(VariablesParsing::press_variable_name(input,physics_name,VariablesParsing::PHYSICS)))
  {
    _press_var.set_is_constraint_var(true);
    this->check_var_subdomain_consistency(_press_var);
  }

  template<unsigned int Dim,typename StrainEnergy>
  void IncompressibleHyperelasticity<Dim,StrainEnergy>::init_context( AssemblyContext & context )
  {
    context.get_element_fe(_press_var.p(),Dim)->get_phi();
    context.get_element_fe(_press_var.p(),Dim)->get_JxW();

    CartesianSolidMechanics<Dim>::init_context(context);
  }

  template<unsigned int Dim,typename StrainEnergy>
  void IncompressibleHyperelasticity<Dim,StrainEnergy>::element_time_derivative
  ( bool compute_jacobian, AssemblyContext & context )
  {
    unsigned int u_var = this->_disp_vars.u();
    unsigned int v_var = this->_disp_vars.v();

    unsigned int w_var = libMesh::invalid_uint;
    if(Dim==3)
      w_var = this->_disp_vars.w();

    unsigned int p_var = this->_press_var.p();

    const MultiphysicsSystem & system = context.get_multiphysics_system();

    // If we have an unsteady solver, then we need to extract the corresponding
    // velocity variable. This allows us to use either a FirstOrderUnsteadySolver
    // or a SecondOrderUnsteadySolver. That is, we get back the velocity variable
    // index for FirstOrderUnsteadySolvers or, if it's a SecondOrderUnsteadySolver,
    // this is actually just giving us back the same variable index.

    // If we only wanted to use a SecondOrderUnsteadySolver, then this
    // step would be unnecessary and we would just
    // populate the _u_var, etc. blocks of the residual and Jacobian.
    unsigned int u_dot_var = system.get_second_order_dot_var(u_var);
    unsigned int v_dot_var = system.get_second_order_dot_var(v_var);
    unsigned int w_dot_var = libMesh::invalid_uint;
    if(Dim==3)
      w_dot_var = system.get_second_order_dot_var(w_var);

    libMesh::DenseSubVector<libMesh::Number> & Fu = context.get_elem_residual(u_dot_var);
    libMesh::DenseSubVector<libMesh::Number> & Fv = context.get_elem_residual(v_dot_var);
    libMesh::DenseSubVector<libMesh::Number> * Fw = nullptr;

    if(Dim==3)
      Fw = &context.get_elem_residual(w_dot_var);

    libMesh::DenseSubMatrix<libMesh::Number> & Kuu = context.get_elem_jacobian(u_dot_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kuv = context.get_elem_jacobian(u_dot_var,v_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvu = context.get_elem_jacobian(v_dot_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number> & Kvv = context.get_elem_jacobian(v_dot_var,v_var);

    libMesh::DenseSubMatrix<libMesh::Number> * Kuw = nullptr;
    libMesh::DenseSubMatrix<libMesh::Number> * Kvw = nullptr;
    libMesh::DenseSubMatrix<libMesh::Number> * Kwu = nullptr;
    libMesh::DenseSubMatrix<libMesh::Number> * Kwv = nullptr;
    libMesh::DenseSubMatrix<libMesh::Number> * Kww = nullptr;

    if(Dim==3)
      {
        Kuw = &context.get_elem_jacobian(u_dot_var,w_var);
        Kvw = &context.get_elem_jacobian(v_dot_var,w_var);
        Kwu = &context.get_elem_jacobian(w_dot_var,u_var);
        Kwv = &context.get_elem_jacobian(w_dot_var,v_var);
        Kww = &context.get_elem_jacobian(w_dot_var,w_var);
      }

    int n_qpoints = context.get_element_qrule().n_points();

    const int n_u_dofs = context.get_dof_indices(u_var).size();

    const std::vector<libMesh::Real> & JxW = this->get_fe(context)->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> > & dphi = this->get_fe(context)->get_dphi();

    IncompressibleHyperelasticityWeakForm<StrainEnergy> weak_form;

    for( int qp=0; qp != n_qpoints; qp++ )
      {
        libMesh::Gradient grad_u, grad_v,grad_w;
        context.interior_gradient(u_var, qp, grad_u);
        context.interior_gradient(v_var, qp, grad_v);

        if(Dim==3)
          context.interior_gradient(w_var, qp, grad_w);

        // After this, all calculations should be dimension independent
        // since the structure of F will be correct for 2D plane strain
        // or 3D
        libMesh::Tensor F;
        if(Dim==2)
          F = this->form_def_gradient(grad_u,grad_v);
        else if(Dim==3)
          F = this->form_def_gradient(grad_u,grad_v,grad_w);



        CartesianHyperlasticity<StrainEnergy> stress_law(F, (*(this->_strain_energy)));

        libMesh::Tensor P(stress_law.pk1_stress());
        const libMesh::Tensor & S = stress_law.get_pk2_stress();

        libMesh::Number J = stress_law.get_J();

        libMesh::Tensor FCinv(F*(stress_law.get_C_inverse()));

        libMesh::Number press;
        context.interior_value(p_var, qp, press);

        for( int i=0; i != n_u_dofs; i++)
          {
            libMesh::RealGradient dphiJ(dphi[i][qp]*JxW[qp]);

            if(Dim==2)
              {
                weak_form.evaluate_internal_stress_residual(P,dphiJ,Fu(i),Fv(i));
                weak_form.evaluate_pressure_stress_residual(J,press,FCinv,dphiJ,Fu(i),Fv(i));
              }
            else if(Dim==3)
              {
                weak_form.evaluate_internal_stress_residual(P,dphiJ,Fu(i),Fv(i),(*Fw)(i));
                weak_form.evaluate_pressure_stress_residual(J,press,FCinv,dphiJ,Fu(i),Fv(i),(*Fw)(i));
              }

            if( compute_jacobian )
              {
                for( int j = 0; j != n_u_dofs; j++ )
                  {
                    if(Dim==2)
                      {
                        weak_form.evaluate_internal_stress_jacobian(S,F,dphiJ,dphi[j][qp],stress_law,
                                                                    Kuu(i,j),Kuv(i,j),
                                                                    Kvu(i,j),Kvv(i,j));
                      }
                    else if(Dim==3)
                      {
                        weak_form.evaluate_internal_stress_jacobian(S,F,dphiJ,dphi[j][qp],stress_law,
                                                                    Kuu(i,j), Kuv(i,j), (*Kuw)(i,j),
                                                                    Kvu(i,j), Kvv(i,j), (*Kvw)(i,j),
                                                                    (*Kwu)(i,j), (*Kwv)(i,j), (*Kww)(i,j));
                      }
                  } // end displacement dof loop

              } // compute jacobian
          }
      }
  }

} // end namespace GRINS
