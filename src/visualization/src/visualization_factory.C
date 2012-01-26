//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010,2011 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
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
