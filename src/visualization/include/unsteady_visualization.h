#ifndef GRINS_UNSTEADY_VISUALIZATION_H
#define GRINS_UNSTEADY_VISUALIZATION_H

#include "visualization.h"

namespace GRINS
{
  class UnsteadyVisualization : public Visualization
  {
  public:

    UnsteadyVisualization(const GetPot& input);
    ~UnsteadyVisualization();

    virtual void output_residual( libMesh::EquationSystems* equation_system,
				  GRINS::MultiphysicsSystem* system,
				  const unsigned int time_step );

  protected:

  };
}
#endif
