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

#ifndef GRINS_ERROR_ESTIMATOR_OPTIONS_H
#define GRINS_ERROR_ESTIMATOR_OPTIONS_H

// C++
#include <string>

// libmMesh forward declarations
class GetPot;

namespace GRINS
{
  //! Container for ErrorEstimator options
  class ErrorEstimatorOptions
  {
  public:

    ErrorEstimatorOptions( const GetPot& input );
    ~ErrorEstimatorOptions(){};

    const std::string& estimator_type() const
    { return _estimator_type; }

    bool patch_reuse() const
    { return _patch_reuse; }

    unsigned char n_adjoint_h_refinements() const
    { return _n_adjoint_h_refinements; }

    unsigned char n_adjoint_p_refinements() const
    { return _n_adjoint_p_refinements; }

    bool compute_qoi_error_estimate() const
    { return _compute_qoi_error_estimate; }

    bool estimator_requires_adjoint() const;

  private:

    void check_dup_input_style( const GetPot& input ) const;

    bool is_old_style( const GetPot& input ) const;

    void parse_old_style(const GetPot& input);

    void parse_new_style(const GetPot& input);

    std::string _estimator_type;
    bool _patch_reuse;
    unsigned char _n_adjoint_h_refinements;
    unsigned char _n_adjoint_p_refinements;
    bool _compute_qoi_error_estimate;

  };

} // end namespace GRINS

#endif // GRINS_ERROR_ESTIMATOR_OPTIONS_H
