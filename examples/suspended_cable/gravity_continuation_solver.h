/*
 * gravity_continuation_solver.h
 *
 *  Created on: Dec 14, 2014
 *      Author: cahaynes
 */

#ifndef GRINS_GRAVITY_CONTINUATION_SOLVER_H_
#define GRINS_GRAVITY_CONTINUATION_SOLVER_H_

//GRINS
#include "grins/grins_steady_solver.h"

namespace GRINS
{
  class GravityContinuationSolver : public SteadySolver
  {
  public:

	  GravityContinuationSolver( const GetPot& input );
    virtual ~GravityContinuationSolver();

    virtual void solve( SolverContext& context );

  protected:

    void increment_gravity( GRINS::MultiphysicsSystem& system,
                             libMesh::Real gravity );

    std::vector<libMesh::Real> _gravity_values;

  };
} // namespace GRINS

#endif /* GRINS_GRAVITY_CONTINUATION_SOLVER_H_ */
