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

// This class
#include "grins/error_estimator_options.h"

// GRINS
#include "grins/common.h"
#include "grins/strategies_parsing.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  ErrorEstimatorOptions::ErrorEstimatorOptions( const GetPot& input )
    : _estimator_type("none"),
      _patch_reuse(false),
      _n_adjoint_h_refinements(1),
      _n_adjoint_p_refinements(0),
      _compute_qoi_error_estimate(false)
  {
    this->check_dup_input_style(input);

    if( this->is_old_style(input) )
      this->parse_old_style(input);
    else
      this->parse_new_style(input);
  }

  void ErrorEstimatorOptions::check_dup_input_style( const GetPot& input ) const
  {
    if( input.have_section("MeshAdaptivity") &&
        input.have_section("Strategies/ErrorEstimation") )
      libmesh_error_msg("ERROR: Cannot use both old and new style of options for ErrorEstimator!");
  }

  bool ErrorEstimatorOptions::is_old_style( const GetPot& input ) const
  {
    return input.have_section("MeshAdaptivity");
  }

  void ErrorEstimatorOptions::parse_old_style(const GetPot& input)
  {
    {
      std::string warning = "WARNING: Using [MeshAdaptivity/<options>] is a DEPRECATED\n";
      warning += "         style of input for ErrorEstimator options. Please\n";
      warning += "         update to use the [Strategies/ErrorEstimation/<options> style.\n";
      grins_warning(warning);
    }

    _estimator_type = input("MeshAdaptivity/estimator_type", "none");
    _patch_reuse = input("MeshAdaptivity/patch_reuse", false);
    _n_adjoint_h_refinements = input("MeshAdaptivity/n_adjoint_h_refinements", 1);
    _n_adjoint_p_refinements = input("MeshAdaptivity/n_adjoint_p_refinements", 0);
    _compute_qoi_error_estimate = input("MeshAdaptivity/compute_qoi_error_estimate", false);
  }

  void ErrorEstimatorOptions::parse_new_style(const GetPot& input)
  {
    _estimator_type = input("Strategies/ErrorEstimation/estimator_type", "none");
    _patch_reuse = input("Strategies/ErrorEstimation/patch_reuse", false);
    _n_adjoint_h_refinements = input("Strategies/ErrorEstimation/n_adjoint_h_refinements", 1);
    _n_adjoint_p_refinements = input("Strategies/ErrorEstimation/n_adjoint_p_refinements", 0);
    _compute_qoi_error_estimate = input("Strategies/ErrorEstimation/compute_qoi_error_estimate", false);
  }

  bool ErrorEstimatorOptions::estimator_requires_adjoint() const
  {
    bool requires_adjoint = false;
    if( _estimator_type == StrategiesParsing::adjoint_residual_error_estimator() ||
        _estimator_type == StrategiesParsing::adjoint_residual_error_estimator() )
      requires_adjoint = true;

    return requires_adjoint;
  }
} // end namespace GRINS
