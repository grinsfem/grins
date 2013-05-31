#include "model_error_assembly.h"
#include "equation_systems.h"

// @TO: for check sum
#include <numeric>

namespace GRINS
{
  // Copied from FEMSystem AssemblyContributions
  // Constructor to set context
  ModelingErrorAssembly::ModelingErrorAssembly(
      MultiphysicsSystem* sys,
      const DofMap* dofindices, // REPL:
      libMesh::ErrorVector& error_vector )
  : _sys(sys),
    _dofindices(dofindices), // REPL:
    _error_vector(error_vector)
    {}

  // Copied from FEMSystem AssemblyContributions::operator()
  // for use with Threads::parallel_for().
  void ModelingErrorAssembly::operator()( const ConstElemRange &range ) const
  {
    // Initialize data structures
    AutoPtr<DiffContext> con = _sys->build_context();
    FEMContext &_fem_context = libmesh_cast_ref<FEMContext&>( *con );
    _sys->init_context( _fem_context );

    // @TO: check residual
    _sys->rhs->zero();

    // The dual solution is assumed unique
    NumericVector<Number> &dual = _sys->get_adjoint_solution( 0 );

    for( ConstElemRange::const_iterator elem_it = range.begin();
         elem_it != range.end();
         ++elem_it )
    {
      const Elem *elem = const_cast<Elem *>( *elem_it );

      // Element dof indices and size
      std::vector<unsigned int> dof_indices, dof_indices_dual;
      // REPL: std::vector<unsigned int> dof_indices;
      _sys->get_dof_map().dof_indices( elem, dof_indices );
      _dofindices->dof_indices( elem, dof_indices_dual ); // REPL:
      unsigned int n_dofs = dof_indices.size();

      // Reinitialization of FEMContext
      _fem_context.pre_fe_reinit( *_sys, elem );
      _fem_context.elem_fe_reinit();

      // Build element-wise residual by calling physics, with 1st argument
      // false as we do not need to compute the jacobian.
      // bool jacobian_computed =
      _sys->time_solver->element_residual( false, _fem_context );

      // Might help
      _sys->get_dof_map().constrain_element_vector(
          _fem_context.elem_residual,
          _fem_context.dof_indices,
          false );

      // Reference to created element residual
      DenseVector<Number> &elem_residual = _fem_context.elem_residual;

      // Reference to element adjoint solution
      DenseVector<Number> elem_dual;

      // @TO: check residual
      _sys->rhs->add_vector( elem_residual, _fem_context.dof_indices );

      // Copied from FEMContext::pre_fe_reinit()
      // Populate element adjoint solution
      elem_dual.resize( n_dofs );
      for( unsigned int i=0; i != n_dofs; ++i )
      {
        elem_dual(i) = dual( dof_indices_dual[i] );
        // REPL: elem_dual(i) = dual( dof_indices[i] );
      }

      // Populate error vector
      unsigned int e_id = (*elem_it)->id();
      for( unsigned int i=0; i != n_dofs; ++i )
      {
        this->_error_vector[e_id] += elem_residual(i) * elem_dual(i);
      }

      // @TO: check sums
      // _sys->rhs->close();
      // if( e_id == 196 )
      // {
      //   elem_residual.print( out );
      //   out << "incorrect sum: " << std::accumulate( _error_vector.begin(),
      //       _error_vector.end(), 0. ) << "\n";
      // }
      // out << "@ elem #" << e_id << " correct sum = " << _sys->rhs->sum() << ", diff in sum = " <<
      //     (_sys->rhs->sum()-std::accumulate( _error_vector.begin(), _error_vector.end(), 0. )) <<
      //     std::endl;
    }
  }
}
