//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
// Copyright (C) 2010-2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

// C++
#include <numeric>
#include <iomanip>

// This class
#include "grins/mesh_adaptive_solver_base.h"
#include "grins/solver_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/mesh_base.h"
#include "libmesh/error_vector.h"

namespace GRINS
{
  MeshAdaptiveSolverBase::MeshAdaptiveSolverBase( const GetPot& input )
    : _error_estimator_options(input),
      _mesh_adaptivity_options(input),
      _refinement_type(INVALID)
  {
    this->set_refinement_type( input, _mesh_adaptivity_options, _refinement_type );
  }

  void MeshAdaptiveSolverBase::build_mesh_refinement( libMesh::MeshBase& mesh )
  {
    _mesh_refinement.reset( new libMesh::MeshRefinement( mesh ) );
    _mesh_refinement->coarsen_by_parents() = _mesh_adaptivity_options.coarsen_by_parents();
    _mesh_refinement->absolute_global_tolerance() = _mesh_adaptivity_options.absolute_global_tolerance();
    _mesh_refinement->nelem_target() = _mesh_adaptivity_options.nelem_target();
    _mesh_refinement->refine_fraction() = _mesh_adaptivity_options.refine_fraction();
    _mesh_refinement->coarsen_fraction() = _mesh_adaptivity_options.coarsen_fraction();
    _mesh_refinement->coarsen_threshold() = _mesh_adaptivity_options.coarsen_threshold();
    _mesh_refinement->node_level_mismatch_limit() = _mesh_adaptivity_options.node_level_mismatch_limit();
    _mesh_refinement->edge_level_mismatch_limit() = _mesh_adaptivity_options.edge_level_mismatch_limit();
    _mesh_refinement->face_level_mismatch_limit() = _mesh_adaptivity_options.face_level_mismatch_limit();
    _mesh_refinement->enforce_mismatch_limit_prior_to_refinement() = _mesh_adaptivity_options.enforce_mismatch_limit_prior_to_refinement();
    _mesh_refinement->max_h_level() = _mesh_adaptivity_options.max_h_level();
  }

  void MeshAdaptiveSolverBase::set_refinement_type( const GetPot& input,
                                                    const MeshAdaptivityOptions& mesh_adaptivity_options,
                                                    MeshAdaptiveSolverBase::RefinementFlaggingType& refinement_type )
  {
    const std::string& refinement_stategy = mesh_adaptivity_options.refinement_strategy();

    if( refinement_stategy == std::string("error_tolerance") )
      {
        refinement_type = ERROR_TOLERANCE;
      }

    else if( refinement_stategy == std::string("nelem_target") )
      {
        refinement_type = N_ELEM_TARGET;

        if( !input.have_variable("MeshAdaptivity/nelem_target") )
          {
            std::cerr << "Error: Must specify nelem_target for" << std::endl
                      << "       for error_tolerance refinement strategy."
                      << std::endl;

            libmesh_error();
          }
      }

    else if( refinement_stategy == std::string("error_fraction") )
      {
        refinement_type = ERROR_FRACTION;
      }

    else if( refinement_stategy == std::string("elem_fraction") )
      {
        refinement_type = ELEM_FRACTION;
      }

    else if( refinement_stategy == std::string("mean_std_dev") )
      {
        refinement_type = MEAN_STD_DEV;
      }

    else
      {
        std::cerr << "Error: Invalid refinement strategy " << refinement_stategy << std::endl
                  << "       Valid refinement strategy options are: absolute_global_tolerance" << std::endl
                  << "                                              error_tolerance" << std::endl
                  << "                                              nelem_target" << std::endl
                  << "                                              error_fraction" << std::endl
                  << "                                              elem_fraction" << std::endl
                  << "                                              mean_std_dev" << std::endl;

        libmesh_error();
      }
  }

  bool MeshAdaptiveSolverBase::check_for_convergence( SolverContext& context,
                                                      const libMesh::ErrorVector& error ) const
  {
    std::cout << "==========================================================" << std::endl
              << "Checking convergence" << std::endl
              << "==========================================================" << std::endl;

    bool converged = false;

    libMesh::Real error_estimate = 0.0;

    if( context.do_adjoint_solve )
      {
        error_estimate = std::accumulate( error.begin(), error.end(), 0.0 );
      }
    else
      {
        error_estimate = error.l2_norm();
      }

    std::cout << "==========================================================" << std::endl
              << "Error estimate = " << error_estimate << std::endl
              << "==========================================================" << std::endl;

    // For now, we just check the norm
    if( std::fabs(error_estimate) <= _mesh_adaptivity_options.absolute_global_tolerance() )
      {
        converged = true;
      }

    return converged;
  }

  void MeshAdaptiveSolverBase::flag_elements_for_refinement( const libMesh::ErrorVector& error )
  {
    switch(_refinement_type)
      {
      case(ERROR_TOLERANCE):
        {
          _mesh_refinement->flag_elements_by_error_tolerance( error );
        }
        break;

      case(N_ELEM_TARGET):
        {
          _mesh_refinement->flag_elements_by_nelem_target( error );
        }
        break;

      case( ERROR_FRACTION ):
        {
          _mesh_refinement->flag_elements_by_error_fraction( error );
        }
        break;

      case( ELEM_FRACTION ):
        {
          _mesh_refinement->flag_elements_by_elem_fraction( error );
        }
        break;

      case( MEAN_STD_DEV ):
        {
          _mesh_refinement->flag_elements_by_mean_stddev( error );
        }
        break;

      case(INVALID):
        {
          std::cerr << "Error: Invalid refinement option!" << std::endl;
          libmesh_error();
        }
        break;

        // Wat?!
      default:
        {
          libmesh_error();
        }

      } // switch(_refinement_type)
  }

  void MeshAdaptiveSolverBase::estimate_error_for_amr( SolverContext& context, libMesh::ErrorVector& error )
  {
    std::cout << "==========================================================" << std::endl
              << "Estimating error" << std::endl
              << "==========================================================" << std::endl;
    context.error_estimator->estimate_error( *context.system, error );

    libMesh::MeshBase& mesh = context.equation_system->get_mesh();

    // Plot error vector
    if( _mesh_adaptivity_options.plot_cell_errors() )
      {
        error.plot_error( _mesh_adaptivity_options.error_plot_prefix()+".exo", mesh );
      }
  }

  void MeshAdaptiveSolverBase::perform_amr( SolverContext& context, const libMesh::ErrorVector& error )
  {
    libMesh::MeshBase& mesh = context.equation_system->get_mesh();

    std::cout << "==========================================================" << std::endl
              << "Performing Mesh Refinement" << std::endl
              << "==========================================================" << std::endl;

    this->flag_elements_for_refinement( error );
    _mesh_refinement->refine_and_coarsen_elements();

    // Dont forget to reinit the system after each adaptive refinement!
    context.equation_system->reinit();

    // This output cannot be toggled in the input file.
    std::cout << "==========================================================" << std::endl
              << "Refined mesh to " << std::setw(12) << mesh.n_active_elem()
              << " active elements" << std::endl
              << "            " << std::setw(16) << context.system->n_active_dofs()
              << " active dofs" << std::endl
              << "==========================================================" << std::endl;
  }

} // end namespace GRINS
