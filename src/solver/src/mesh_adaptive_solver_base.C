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

// This class
#include "grins/mesh_adaptive_solver_base.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/mesh_base.h"

namespace GRINS
{
  MeshAdaptiveSolverBase::MeshAdaptiveSolverBase( const GetPot& input )
    : Solver( input ),
      _max_r_steps( input("MeshAdaptivity/max_r_steps", 5) ),
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
      _refinement_type(INVALID)
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
    // Check that either nelem_target or global tolerance was set
    if( !input.have_variable("MeshAdaptivity/absolute_global_tolerance") &&
        !input.have_variable("MeshAdaptivity/nelem_target") )
      {
        std::cerr << "Error: Must specify either global tolerance or element number target" << std::endl
                  << "       for mesh adaptive solver." << std::endl;
        libmesh_error();
      }

    // Make sure *both* weren't set
    if( input.have_variable("MeshAdaptivity/absolute_global_tolerance") &&
        input.have_variable("MeshAdaptivity/nelem_target") )
      {
        std::cerr << "Error: Can only specify either global tolerance or element number target" << std::endl
                  << "       for mesh adaptive solver." << std::endl;
        libmesh_error();
      }

    if( input.have_variable("MeshAdaptivity/absolute_global_tolerance") )
      {
        refinement_type = ERROR_TOLERANCE;
      }
    
    if( input.have_variable("MeshAdaptivity/nelem_target") )
      {
        refinement_type = N_ELEM_TARGET;
      }

    return;
  }

  bool MeshAdaptiveSolverBase::check_for_adjoint_solve( const GetPot& input )
  {
    std::string error_estimator = input("MeshAdaptivity/estimator_type");

    bool do_adjoint_solve = false;

    if( error_estimator.find("adjoint") != error_estimator.end() )
      {
        do_adjoint_solve = true;
      }

    return do_adjoint_solve;
  }

  bool MeshAdaptiveSolverBase::check_for_convergence()
  {
    bool converged = false;

    switch(_refinement_type)
      {
      case(ERROR_TOLERANCE):
        {

        }
        break;

      case(N_ELEM_TARGET):
        {
        }
        break;
        
      case(INVALID):
        {

        }
        break;

      // Wat?!
      default:
        {
          libmesh_error();
        }

      } // switch(_refinement_type)
    

    if( this->_absolute_global_tolerance >= 0. && this->_nelem_target == 0.)
          {
            _mesh_refinement->flag_elements_by_error_tolerance( error );
          }
        // Adaptively refine based on reaching a target number of elements
        else
          {
            if( mesh.n_active_elem() >= this->_nelem_target )
              {
                std::cout<<"We reached the target number of elements."<<std::endl <<std::endl;
                break;
              }
            
            _mesh_refinement->flag_elements_by_nelem_target( error );
          }

    return converged;
  }

} // end namespace GRINS
