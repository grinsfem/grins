#include "model_error_estimator.h"

namespace GRINS
{
  ModelErrorEstimator::ModelErrorEstimator(
    GRINS::MultiphysicsSystem* system,
    libMesh::MeshBase* mesh,
    libMesh::ErrorVector* error_vector )
  : _system( system ),
    _mesh( mesh ),
    _error_vector( *error_vector )
    {}

  ModelErrorEstimator::~ModelErrorEstimator(){}

  void ModelErrorEstimator::assembly()
  {
    // Copied from FEMSystem::assembly()
    // Update ghost cells, their data, etc
    _system->update();
  
    // Reset the residual vector
    _system->rhs->zero();

    // Reset error_vector (just in case mesh changed)
    this->_error_vector.clear();
    const unsigned int n_elems = _mesh->n_elem();
    this->_error_vector.resize( n_elems );
  
    // Parallel assemble element residual, gets passed:
    // processor element range
    // functor object ModelingErrorAssembly
    //ElemRange elem_range;
    //Threads::parallel_for( 
    //  elem_range.reset(
    //  _mesh->active_local_elements_begin(),
    //  _mesh->active_local_elements_end() );
    //  ModelingErrorAssembly assembly( _system, _error_vector ) );

    libMesh::MeshBase::element_iterator iter_begin = _mesh->active_local_elements_begin();
    libMesh::MeshBase::element_iterator iter_end = _mesh->active_local_elements_end();
    ModelingErrorAssembly assembly( _system, _error_vector );
    assembly( iter_begin, iter_end );
  
    _system->rhs->close();

    return;
  }
}
