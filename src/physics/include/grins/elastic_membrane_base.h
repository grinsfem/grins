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

#ifndef GRINS_ELASTIC_MEMBRANE_BASE_H
#define GRINS_ELASTIC_MEMBRANE_BASE_H

//GRINS
#include "grins/curvilinear_solid_mechanics.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  class ElasticMembraneBase : public CurvilinearSolidMechanics<2>
  {
  public:

    ElasticMembraneBase( const PhysicsName& physics_name,
                         const GetPot& input,
                         bool is_compressible);

    virtual ~ElasticMembraneBase() = default;

    virtual void init_variables( libMesh::FEMSystem* system ) override;

  protected:

    //! Implementation of mass_residual.
    /*! The mu argument is needed to support Rayleigh Damping. */
    void mass_residual_impl( bool compute_jacobian,
                             AssemblyContext& context,
                             InteriorFuncType interior_solution,
                             VarDerivType get_solution_deriv,
                             libMesh::Real mu = 1.0 );

    void compute_metric_tensors( unsigned int qp,
                                 const libMesh::FEBase& elem,
                                 const AssemblyContext& context,
                                 const libMesh::Gradient& grad_u,
                                 const libMesh::Gradient& grad_v,
                                 const libMesh::Gradient& grad_w,
                                 libMesh::TensorValue<libMesh::Real>& a_cov,
                                 libMesh::TensorValue<libMesh::Real>& a_contra,
                                 libMesh::TensorValue<libMesh::Real>& A_cov,
                                 libMesh::TensorValue<libMesh::Real>& A_contra,
                                 libMesh::Real& lambda_sq);

    StressStrainLaw _stress_strain_law;

    bool _is_compressible;

    //! Membrane thickness
    libMesh::Real _h0;

    //! Variable index for lambda_sq variable
    VariableIndex _lambda_sq_var;

  };

}

#endif // GRINS_ELASTIC_MEMBRANE_BASE_H
