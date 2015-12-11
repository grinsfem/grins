//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_ELASTIC_CABLE_H
#define GRINS_ELASTIC_CABLE_H


//GRINS
#include "grins/elastic_cable_base.h"

//LIBMESH
#include "libmesh/fe_base.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  class ElasticCable : public ElasticCableBase
  {
  public:

    ElasticCable( const PhysicsName& physics_name, const GetPot& input,
                  bool lambda_sq_var );

    virtual ~ElasticCable();

    //! Register postprocessing variables for ElasticCable
    virtual void register_postprocessing_vars( const GetPot& input,
                                               PostProcessedQuantities<libMesh::Real>& postprocessing );

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context,
                                          CachedValues& cache );

    virtual void side_time_derivative( bool compute_jacobian,
                                       AssemblyContext& context,
                                       CachedValues& cache );

    virtual void mass_residual( bool compute_jacobian,
                                AssemblyContext& context,
                                CachedValues& cache );

    //! Compute the registered postprocessed quantities
    virtual void compute_postprocessed_quantity( unsigned int quantity_index,
                                                 const AssemblyContext& context,
                                                 const libMesh::Point& point,
                                                 libMesh::Real& value );

  private:

    ElasticCable();

    // This is straight up copied from libMesh. Should make this a friend or public.
    libMesh::AutoPtr<libMesh::FEGenericBase<libMesh::Real> > build_new_fe( const libMesh::Elem& elem,
                                                                           const libMesh::FEGenericBase<libMesh::Real>* fe,
                                                                           const libMesh::Point p );

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

    //! Index from registering this quantity. Each component will have it's own index.
    std::vector<unsigned int> _stress_indices;

    //! Index from registering this quantity. Each component will have it's own index.
    std::vector<unsigned int> _strain_indices;

    //! Index from registering this quantity. Each component will have it's own index.
    std::vector<unsigned int> _force_indices;

  };

} // end namespace GRINS


#endif /* GRINS_ELASTIC_CABLE_H */
