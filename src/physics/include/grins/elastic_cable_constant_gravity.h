/*
 * elastic_cable_constant_gravity.h
 *
 *  Created on: Dec 14, 2014
 *      Author: cahaynes
 */

#ifndef GRINS_ELASTIC_CABLE_CONSTANT_GRAVITY_H_
#define GRINS_ELASTIC_CABLE_CONSTANT_GRAVITY_H_

//GRINS
#include "grins/physics.h"
#include "grins/solid_mechanics_fe_variables.h"

namespace GRINS
{
  class ElasticCableConstantGravity : public Physics
  {
  public:
	  ElasticCableConstantGravity( const GRINS::PhysicsName& physics_name, const GetPot& input );
    virtual ~ElasticCableConstantGravity();

    //! Initialize variables for this physics.
    virtual void init_variables( libMesh::FEMSystem* system );

    virtual void set_time_evolving_vars( libMesh::FEMSystem* system );

    //! Initialize context for added physics variables
    virtual void init_context( AssemblyContext& context );

    //! Time dependent part(s) of physics for element interiors
    virtual void element_time_derivative( bool compute_jacobian,
                                          AssemblyContext& context,
                                          CachedValues& cache );

    void reset_gravity( libMesh::Real gravity_in );

  protected:

    SolidMechanicsFEVariables _disp_vars;

  private:

    ElasticCableConstantGravity();

    libMesh::Real _gravity;
    libMesh::Real _A0;
    libMesh::Real _rho0;
  };

} // end namespace GRINS

#endif /* GRINS_ELASTIC_CABLE_CONSTANT_GRAVITY_H_ */
