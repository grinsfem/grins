//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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

// This class
#include "grins/mesh_adaptive_solver_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/mesh_base.h"
#include "libmesh/error_vector.h"

namespace GRINS
{
  MeshAdaptiveSolverBase::MeshAdaptiveSolverBase( const GetPot& input )
    : Solver( input ),
      _max_refinement_steps( input("MeshAdaptivity/max_refinement_steps", 5) ),
      _coarsen_by_parents(true),
      _absolute_global_tolerance( input("MeshAdaptivity/absolute_global_tolerance", 0) ),
      _nelem_target( input("MeshAdaptivity/nelem_target", 0) ),
      _refine_fraction( input("MeshAdaptivity/refine_percentage", 0.2) ),
      _coarsen_fraction( input("MeshAdaptivity/coarsen_percentage", 0.2) ),
      _coarsen_threshold( input("MeshAdaptivity/coarsen_threshold", 0) ),
      _output_adjoint_sol( input("MeshAdaptivity/output_adjoint_sol", false) ),
      _plot_cell_errors( input("MeshAdaptivity/plot_cell_errors", false) ),
      _error_plot_prefix( input("MeshAdaptivity/error_plot_prefix", "cell_error") ),
      _do_adjoint_solve(false),
      _refinement_type(INVALID),
      _mesh_refinement(NULL)
  {
    this->set_refinement_type( input, _refinement_type );

    _do_adjoint_solve = this->check_for_adjoint_solve( input );

    return;
  }
  
  MeshAdaptiveSolverBase::~MeshAdaptiveSolverBase()
  {
    return;
  }

  void MeshAdaptiveSolverBase::build_mesh_refinement( libMesh::MeshBase& mesh )
  {
    _mesh_refinement.reset( new libMesh::MeshRefinement( mesh ) );
    _mesh_refinement->coarsen_by_parents() = _coarsen_by_parents;
    _mesh_refinement->absolute_global_tolerance() = _absolute_global_tolerance;
    _mesh_refinement->nelem_target() = _nelem_target;
    _mesh_refinement->refine_fraction() = _refine_fraction;
    _mesh_refinement->coarsen_fraction() = _coarsen_fraction;  
    _mesh_refinement->coarsen_threshold() = _coarsen_threshold;

    return; 
  }

  void MeshAdaptiveSolverBase::set_refinement_type( const GetPot& input,
                                                    MeshAdaptiveSolverBase::RefinementFlaggingType& refinement_type )
  {
    std::string refinement_stategy = input("MeshAdaptivity/refinement_strategy", "elem_fraction" ); 

    if( !input.have_variable("MeshAdaptivity/absolute_global_tolerance" ) )
      {
        std::cerr << "Error: Must specify absolute_global_tolerance for" << std::endl
                  << "       adaptive refinement algorithms."
                  << std::endl;
        
        libmesh_error();
      }

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

    return;
  }

  bool MeshAdaptiveSolverBase::check_for_adjoint_solve( const GetPot& input ) const
  {
    std::string error_estimator = input("MeshAdaptivity/estimator_type", "none");

    bool do_adjoint_solve = false;

    if( error_estimator.find("adjoint") != std::string::npos )
      {
        do_adjoint_solve = true;
      }

    return do_adjoint_solve;
  }

  bool MeshAdaptiveSolverBase::check_for_convergence( const libMesh::ErrorVector& error ) const
  {
    bool converged = false;

    libMesh::Real error_estimate = 0.0;

    if( _do_adjoint_solve )
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
    if( error_estimate <= _absolute_global_tolerance )
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

    return;
  }

} // end namespace GRINS
