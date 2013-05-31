#ifndef MODEL_ERROR_ASSEMBLY_H
#define MODEL_ERROR_ASSEMBLY_H

// GRINS
#include "multiphysics_sys.h"
#include "solver_factory.h"

// libMesh
#include "dof_map.h"
#include "diff_system.h"
#include "elem.h"
#include "fe_base.h"
#include "fem_context.h"
#include "fem_system.h"
#include "mesh_base.h"
#include "numeric_vector.h"
#include "elem_range.h"

namespace GRINS
{
  // Copied form FEMSystem AssemblyContributions
  class ModelingErrorAssembly
  {
  public:

    // Constructor to set context
    ModelingErrorAssembly(
        MultiphysicsSystem* sys,
        const DofMap *dofindices, // REPL:
        libMesh::ErrorVector& error_vector );
  
    // operator() for use with Threads::parallel_for().
    void operator()( const ConstElemRange &range ) const;
  
  private:

    // This needs to be a standard pointer, as _equation_system will own and destroy the object.
    MultiphysicsSystem* _sys;
    const DofMap *_dofindices; // REPL:
    libMesh::ErrorVector&  _error_vector;
  };
}
#endif // MODEL_ERROR_ASSEMBLY_H
