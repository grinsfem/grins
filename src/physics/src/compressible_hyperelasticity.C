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
#include "grins/compressible_hyperelasticity.h"

// GRINS
#include "grins/materials_parsing.h"
#include "grins/multiphysics_sys.h"

// libMesh
#include "libmesh/quadrature.h"

namespace GRINS
{
  template<typename StrainEnergy>
  CompressibleHyperelasticity<StrainEnergy>::CompressibleHyperelasticity
  ( const PhysicsName & physics_name, const GetPot & input )
    : ThreeDSolidMechanicsBase(physics_name,PhysicsNaming::compressible_hyperelasticity(),input),
      _strain_energy(nullptr)
  {
    const std::string material =
      MaterialsParsing::material_name(input,PhysicsNaming::compressible_hyperelasticity());

    _strain_energy.reset(new StrainEnergy(input,material));
  }

  template<typename StrainEnergy>
  void CompressibleHyperelasticity<StrainEnergy>::element_time_derivative( bool compute_jacobian,
                                                                           AssemblyContext & context )
  {
    unsigned int u_var = this->_disp_vars.u();
    unsigned int v_var = this->_disp_vars.v();
    unsigned int w_var = this->_disp_vars.w();

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
    unsigned int w_dot_var = system.get_second_order_dot_var(w_var);

    libMesh::DenseSubVector<libMesh::Number> & Fu = context.get_elem_residual(u_dot_var);
    libMesh::DenseSubVector<libMesh::Number> & Fv = context.get_elem_residual(v_dot_var);
    libMesh::DenseSubVector<libMesh::Number> & Fw = context.get_elem_residual(w_dot_var);

    libMesh::DenseSubMatrix<libMesh::Number>& Kuu = context.get_elem_jacobian(u_dot_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number>& Kuv = context.get_elem_jacobian(u_dot_var,v_var);
    libMesh::DenseSubMatrix<libMesh::Number>& Kuw = context.get_elem_jacobian(u_dot_var,w_var);

    libMesh::DenseSubMatrix<libMesh::Number>& Kvu = context.get_elem_jacobian(v_dot_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number>& Kvv = context.get_elem_jacobian(v_dot_var,v_var);
    libMesh::DenseSubMatrix<libMesh::Number>& Kvw = context.get_elem_jacobian(v_dot_var,w_var);

    libMesh::DenseSubMatrix<libMesh::Number>& Kwu = context.get_elem_jacobian(w_dot_var,u_var);
    libMesh::DenseSubMatrix<libMesh::Number>& Kwv = context.get_elem_jacobian(w_dot_var,v_var);
    libMesh::DenseSubMatrix<libMesh::Number>& Kww = context.get_elem_jacobian(w_dot_var,w_var);

    int n_qpoints = context.get_element_qrule().n_points();

    const int n_u_dofs = context.get_dof_indices(u_var).size();

    const std::vector<libMesh::Real> & JxW = this->get_fe(context)->get_JxW();

    const std::vector<std::vector<libMesh::RealGradient> > & dphi = this->get_fe(context)->get_dphi();

    const int dim = 3;

    for( int qp=0; qp != n_qpoints; qp++ )
      {
        libMesh::Gradient grad_u, grad_v,grad_w;
        context.interior_gradient(u_var, qp, grad_u);
        context.interior_gradient(v_var, qp, grad_v);
        context.interior_gradient(w_var, qp, grad_w);

        libMesh::Tensor F = this->form_def_gradient(grad_u,grad_v,grad_w);
        libMesh::Tensor C = this->compute_right_cauchy_def(F);

        libMesh::Number I1, I2, I3;
        this->compute_invariants(C,I1,I2,I3);

        libMesh::Number dWdI1 = _strain_energy->dI1(I1,I2,I3);
        libMesh::Number dWdI2 = _strain_energy->dI2(I1,I2,I3);
        libMesh::Number dWdI3 = _strain_energy->dI3(I1,I2,I3);

        libMesh::Tensor Cinv = C.inverse();

        libMesh::Tensor S = this->compute_pk2_stress(C,Cinv,I1,I3,dWdI1,dWdI2,dWdI3);

        libMesh::Tensor P(F*S);

        for( int i=0; i != n_u_dofs; i++)
          {
            for( int alpha = 0; alpha < dim; alpha++)
              {
                libMesh::Real dphiJ = dphi[i][qp](alpha)*JxW[qp];

                Fu(i) += P(0,alpha)*dphiJ;
                Fv(i) += P(1,alpha)*dphiJ;
                Fw(i) += P(2,alpha)*dphiJ;
              }

            if( compute_jacobian )
              {
                libmesh_not_implemented();
              }
          }
      }
  }

  template<typename StrainEnergy>
  libMesh::Tensor CompressibleHyperelasticity<StrainEnergy>::compute_pk2_stress
  ( const libMesh::Tensor & C,
    const libMesh::Tensor & Cinv,
    const libMesh::Number & I1,
    const libMesh::Number & I3,
    const libMesh::Number & dWdI1,
    const libMesh::Number & dWdI2,
    const libMesh::Number & dWdI3 ) const
  {
    const int dim = 3;
    libMesh::Tensor S;

    for( int i = 0; i < dim; i++ )
      {
        for (int j = 0; j < dim; j++ )
          {
            libMesh::Real dij = this->delta(i,j);
            S(i,j) = 2*( dWdI1*dij + dWdI2*(I1*dij-C(i,j)) + dWdI3*(I3*Cinv(i,j)) );
          }
      }

    return S;
  }
} // end namespace GRINS
