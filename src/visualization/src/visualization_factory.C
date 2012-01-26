//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "visualization_factory.h"

GRINS::VisualizationFactory::VisualizationFactory( const GetPot& input )
  : _transient( input("unsteady-solver/transient", false ) ),
    _input( input )
{
  return;
}

GRINS::VisualizationFactory::~VisualizationFactory()
{
  return;
}

libmesh::AutoPtr<GRINS::Visualization> 
GRINS::VisualizationFactory::build()
{
  GRINS::Visualization* vis;

  if(_transient) 
    vis = new UnsteadyVisualization( _input );
  else
    vis = new SteadyVisualization( _input );

  return libMesh::AutoPtr<GRINS::Visualization>( vis );
}
