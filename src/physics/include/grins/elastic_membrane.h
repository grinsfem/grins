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

#ifndef GRINS_ELASTIC_MEMBRANE_H
#define GRINS_ELASTIC_MEMBRANE_H

//GRINS
#include "grins/elastic_membrane_base.h"
#include "grins/elasticity_tensor.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  class ElasticMembrane : public ElasticMembraneBase<StressStrainLaw>
  {
  public:

    ElasticMembrane( const GRINS::PhysicsName& physics_name, const GetPot& input,
                     bool is_compressible );

    virtual ~ElasticMembrane(){};

    //! Register postprocessing variables for ElasticMembrane
    virtual void register_postprocessing_vars( const GetPot& input,
                                               PostProcessedQuantities<libMesh::Real> & postprocessing );

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext & context );

    virtual void element_constraint( bool compute_jacobian,
                                     AssemblyContext & context );

    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext & context )
    { this->mass_residual_impl(compute_jacobian,
                               context,
                               &libMesh::FEMContext::interior_accel,
                               &libMesh::DiffContext::get_elem_solution_accel_derivative); }

    //! Compute the registered postprocessed quantities
    virtual void compute_postprocessed_quantity( unsigned int quantity_index,
                                                 const AssemblyContext& context,
                                                 const libMesh::Point& point,
                                                 libMesh::Real& value );

    //! Precompute data needed for get_stress inline function
    void get_grad_disp( const AssemblyContext & context,
                        unsigned int qp,
                        libMesh::Gradient & grad_u,
                        libMesh::Gradient & grad_v,
                        libMesh::Gradient & grad_w );


    //! Precompute stress and elasticity
    void get_stress_and_elasticity( const AssemblyContext & context,
                                    unsigned int qp,
                                    const libMesh::Gradient & grad_u,
                                    const libMesh::Gradient & grad_v,
                                    const libMesh::Gradient & grad_w,
                                    libMesh::TensorValue<libMesh::Real> & tau,
                                    ElasticityTensor & C );

  private:

    ElasticMembrane();

    //! Index from registering this quantity for postprocessing. Each component will have it's own index.
    std::vector<unsigned int> _stress_indices;

    //! Index from registering sigma_zz for postprocessing. Mainly for sanity checking.
    unsigned int _stress_zz_index;

    //! Index from registering this quantity for postprocessing. Each component will have it's own index.
    std::vector<unsigned int> _strain_indices;

  }; //end class ElasticMembrane

} // end namespace GRINS

#endif // GRINS_ELASTIC_MEMBRANE_H
