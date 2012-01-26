#ifndef GRINS_UNSTEADY_VISUALIZATION_H
#define GRINS_UNSTEADY_VISUALIZATION_H

// libMesh
#include "auto_ptr.h"

// GRINS
#include "visualization.h"

namespace GRINS
{
  class UnsteadyVisualization : public Visualization
  {
  public:

    UnsteadyVisualization(const GetPot& input);
    ~UnsteadyVisualization();

    virtual void output_residual( std::tr1::shared_ptr<libMesh::EquationSystems> equation_system,
				  GRINS::MultiphysicsSystem* system,
				  const unsigned int time_step,
				  const Real time);

  protected:

  };
}
#endif
