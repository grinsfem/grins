//-----------------------------------------------------------------------bl-
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#ifndef SOLVER_FACTORY_H
#define SOLVER_FACTORY_H

//libMesh
#include "auto_ptr.h"
#include "getpot.h"

//GRINS
#include "grins_solver.h"
#include "steady_solver.h"
#include "unsteady_solver.h"

namespace GRINS
{
  //! This object handles constructing the solver to be used.
  /*! To allow the user to easily extend the (limited) available solvers,
      the solver construction is handled in this object. Note that a
      libMesh::AutoPtr is returned to transfer ownership away from this class.
  */
  class SolverFactory
  {
  public:

    SolverFactory(const GetPot& input);
    virtual ~SolverFactory();

    virtual void read_input_options( const GetPot& input );

    //! Builds GRINS::Solver object.
    /*! Users should override this method to construct 
        their own solvers. Note that a libMesh::AutoPtr is 
	returned to transfer ownership away from this class. */
    virtual libMesh::AutoPtr<GRINS::Solver> build();

  protected:

    //! All we need to distinguish between steady and unsteady solver.
    bool _transient;
  };
}
#endif //SOLVER_FACTORY_H
