//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_VISUALIZATION_FACTORY_H
#define GRINS_VISUALIZATION_FACTORY_H

#include "boost/tr1/memory.hpp"

// libMesh
#include "libmesh/getpot.h"

// GRINS
#include "grins/visualization.h"
#include "grins/steady_visualization.h"
#include "grins/unsteady_visualization.h"

namespace GRINS
{
  class VisualizationFactory
  {
  public:

    VisualizationFactory();
    ~VisualizationFactory();

    virtual std::tr1::shared_ptr<GRINS::Visualization> build(const GetPot& input);
  };
}
#endif //GRINS_VISUALIZATION_FACTORY_H
