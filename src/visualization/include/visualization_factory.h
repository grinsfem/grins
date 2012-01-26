//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef GRINS_VISUALIZATION_FACTORY_H
#define GRINS_VISUALIZATION_FACTORY_H

// libMesh
#include "getpot.h"
#include "auto_ptr.h"

// GRINS
#include "visualization.h"
#include "steady_visualization.h"
#include "unsteady_visualization.h"

namespace GRINS
{
  class VisualizationFactory
  {
  public:

    VisualizationFactory( const GetPot& input );
    ~VisualizationFactory();

    virtual libmesh::AutoPtr<GRINS::Visualization> build();

  protected:

    bool _transient;
    
    const GetPot& _input;
  };
}
#endif //GRINS_VISUALIZATION_FACTORY_H
