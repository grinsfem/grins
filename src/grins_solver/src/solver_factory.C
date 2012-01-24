//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "solver_factory.h"

GRINS::SolverFactory::SolverFactory(const GetPot& input)
  :_transient(false)
{
  this->read_input_options(input);
  return;
}

GRINS::SolverFactory::~SolverFactory()
{
  return;
}

void GRINS::SolverFactory::read_input_options( const GetPot& input )
{
  _transient = input("unsteady-solver/transient", false );
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
