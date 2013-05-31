#ifndef MODEL_ERROR_ESTIMATOR_H
#define MODEL_ERROR_ESTIMATOR_H

// libMesh
#include "dof_map.h"
#include "elem.h"
#include "equation_systems.h"
#include "fe_base.h"
#include "fem_context.h"
#include "fem_system.h"
#include "libmesh_config.h"
#include "libmesh_logging.h"
#include "location_maps.h"
#include "mesh_base.h"
#include "numeric_vector.h"
#include "parallel.h"
#include "point_locator_base.h"
#include "quadrature.h"

// GRINS
#include "multiphysics_sys.h"
#include "model_error_assembly.h"

namespace GRINS
{
  class ModelErrorEstimator
  {
    public:

      ModelErrorEstimator(
        MultiphysicsSystem* system,
        MeshBase* mesh,
        ErrorVector* error_vector );

      ~ModelErrorEstimator();

      void assembly();

    private:
      // This needs to be a standard pointer, as _equation_system will own and destroy the object.
      MultiphysicsSystem* _system;
      MeshBase* _mesh;
      libMesh::ErrorVector& _error_vector;
  };
} // namespace GRINS
#endif // MODEL_ERROR_ESTIMATOR_H
