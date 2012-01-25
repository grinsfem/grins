//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "solver_factory.h"

GRINS::SolverFactory::SolverFactory(const GetPot& input)
  :_transient( input("unsteady-solver/transient", false ) )
{
  return;
}

GRINS::SolverFactory::~SolverFactory()
{
  return;
}

libMesh::AutoPtr<GRINS::Solver> build()
{
  GRINS::Solver* solver;
  if(_transient)
    {
      solver = new GRINS::UnsteadySolver();
    }
  else
    {
      solver = new GRINS::SteadySolver();
    }
  return libMesh::AutoPtr(solver);
}
