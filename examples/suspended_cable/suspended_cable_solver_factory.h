/*
 * suspended_cable_solver_factory.h
 *
 *  Created on: Dec 14, 2014
 *      Author: cahaynes
 */

#ifndef GRINS_SUSPENDED_CABLE_SOLVER_FACTORY_H_
#define GRINS_SUSPENDED_CABLE_SOLVER_FACTORY_H_

//GRINS
#include "grins/solver_factory.h"

namespace GRINS
{
  class SuspendedCableSolverFactory : public SolverFactory
  {
  public:

	SuspendedCableSolverFactory(){};
    virtual ~SuspendedCableSolverFactory(){};

    virtual std::tr1::shared_ptr<GRINS::Solver> build(const GetPot& input);
  };

} // end namespace GRINS

#endif /* GRINS_SUSPENDED_CABLE_SOLVER_FACTORY_H_ */
