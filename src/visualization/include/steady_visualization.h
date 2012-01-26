//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_STEADY_VISUALIZATION_H
#define GRINS_STEADY_VISUALIZATION_H

#include "visualization.h"

namespace GRINS
{
  class SteadyVisualization : public Visualization
  {
  public:

    SteadyVisualization(const GetPot& input);
    ~SteadyVisualization();

    virtual void output_residual( libMesh::EquationSystems* equation_system,
				  GRINS::MultiphysicsSystem* system,
				  const unsigned int time_step );

  protected:

  };
}
#endif
