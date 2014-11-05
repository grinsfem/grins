//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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
#include "grins/physics.h"
#include "grins/solid_mechanics_fe_variables.h"

namespace GRINS
{
  template<typename StressStrainLaw>
  class ElasticMembrane : public Physics
  {
  public:

    ElasticMembrane( const GRINS::PhysicsName& physics_name, const GetPot& input,
                     bool lambda_sq_coupled, bool lambda_sq_var );
    virtual ~ElasticMembrane();

    //! Initialize variables for this physics.
    virtual void init_variables( libMesh::FEMSystem* system );

    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

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
    
  protected:

    SolidMechanicsFEVariables _disp_vars;


  private:

    ElasticMembrane();

    StressStrainLaw _stress_strain_law;

    bool _lambda_sq_coupled;
    bool _lambda_sq_var;

  };

} // end namespace GRINS

#endif // GRINS_ELASTIC_MEMBRANE_H
