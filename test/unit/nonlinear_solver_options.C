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

#include "grins_config.h"

#ifdef GRINS_HAVE_CPPUNIT

#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

#include "grins/nonlinear_solver_options.h"

namespace GRINSTesting
{
  class NonlinearSolverOptionsTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( NonlinearSolverOptionsTest );

    CPPUNIT_TEST( test_defaults );
    CPPUNIT_TEST( test_eval );

    CPPUNIT_TEST_SUITE_END();

  protected:

    GetPot _input;

  public:

    void setup_nondefault_inputfile()
    {
      const char* text =
        "[linear-nonlinear-solver]                           \n"
        "max_nonlinear_iterations = '1'                      \n"
        "relative_step_tolerance = '1.0e-9'                  \n"
        "absolute_step_tolerance = '1.0e-10'                 \n"
        "relative_residual_tolerance = '1.0e-10'             \n"
        "absolute_residual_tolerance = '1.0e-13'             \n"
        "continue_after_backtrack_failure  = 'true'          \n"
        "continue_after_max_iterations  = 'true'             \n"
        "require_residual_reduction  = 'false'               \n"
        "verify_analytic_jacobians  = '1.0e-6'               \n"
        "use_numerical_jacobians_only  = 'true'              \n"
        "numerical_jacobian_h  = '1.0e-3'                    \n"
        "numerical_jacobian_h_variables = 'U y Vz'           \n"
        "numerical_jacobian_h_values = '1.0  2.0e-6  3.0e-8' \n"
        "type = 'libmesh_petsc_diff'                         \n"
        "[]\n";

      std::stringstream inputfile;
      inputfile << text;

      _input.parse_input_stream(inputfile);
    }

    //! Capture any changes to the default values
    void test_defaults()
    {
      GRINS::NonlinearSolverOptions options(_input);

      CPPUNIT_ASSERT_EQUAL( 10, (int)options.max_nonlinear_iterations() );
      CPPUNIT_ASSERT_EQUAL( 1.0e-6, options.relative_step_tolerance() );
      CPPUNIT_ASSERT_EQUAL( 0.0, options.absolute_step_tolerance() );
      CPPUNIT_ASSERT_EQUAL( 1.0e-15, options.relative_residual_tolerance() );
      CPPUNIT_ASSERT_EQUAL( 0.0, options.absolute_residual_tolerance() );
      CPPUNIT_ASSERT_EQUAL( false, options.continue_after_backtrack_failure() );
      CPPUNIT_ASSERT_EQUAL( false, options.continue_after_max_iterations() );
      CPPUNIT_ASSERT_EQUAL( true, options.require_residual_reduction() );
      CPPUNIT_ASSERT_EQUAL( 0.0, options.verify_analytic_jacobians() );
      CPPUNIT_ASSERT_EQUAL( false, options.use_numerical_jacobians_only() );
      CPPUNIT_ASSERT_EQUAL( libMesh::TOLERANCE, options.numerical_jacobian_h() );

      std::vector<std::string> variables;
      std::vector<libMesh::Real> values;
      options.numerical_jacobian_h_vars_and_vals(variables,values);
      CPPUNIT_ASSERT(variables.empty());
      CPPUNIT_ASSERT(values.empty());

      CPPUNIT_ASSERT_EQUAL(GRINS::DiffSolverNames::newton_solver(), options.type() );
    }

    void test_eval()
    {
      this->setup_nondefault_inputfile();
      GRINS::NonlinearSolverOptions options(_input);

      CPPUNIT_ASSERT_EQUAL( 1, (int)options.max_nonlinear_iterations() );
      CPPUNIT_ASSERT_EQUAL( 1.0e-9, options.relative_step_tolerance() );
      CPPUNIT_ASSERT_EQUAL( 1.0e-10, options.absolute_step_tolerance() );
      CPPUNIT_ASSERT_EQUAL( 1.0e-10, options.relative_residual_tolerance() );
      CPPUNIT_ASSERT_EQUAL( 1.0e-13, options.absolute_residual_tolerance() );
      CPPUNIT_ASSERT_EQUAL( true, options.continue_after_backtrack_failure() );
      CPPUNIT_ASSERT_EQUAL( true, options.continue_after_max_iterations() );
      CPPUNIT_ASSERT_EQUAL( false, options.require_residual_reduction() );
      CPPUNIT_ASSERT_EQUAL( 1.e-6, options.verify_analytic_jacobians() );
      CPPUNIT_ASSERT_EQUAL( true, options.use_numerical_jacobians_only() );
      CPPUNIT_ASSERT_EQUAL( 1.0e-3, options.numerical_jacobian_h() );

      std::vector<std::string> variables;
      std::vector<libMesh::Real> values;
      options.numerical_jacobian_h_vars_and_vals(variables,values);
      CPPUNIT_ASSERT_EQUAL(3, (int)variables.size());
      CPPUNIT_ASSERT_EQUAL(std::string("U"),variables[0]);
      CPPUNIT_ASSERT_EQUAL(std::string("y"),variables[1]);
      CPPUNIT_ASSERT_EQUAL(std::string("Vz"),variables[2]);

      CPPUNIT_ASSERT_EQUAL(3, (int)values.size());
      CPPUNIT_ASSERT_EQUAL(1.0,values[0]);
      CPPUNIT_ASSERT_EQUAL(2.0e-6,values[1]);
      CPPUNIT_ASSERT_EQUAL(3.0e-8,values[2]);

      CPPUNIT_ASSERT_EQUAL(GRINS::DiffSolverNames::petsc_diff_solver(), options.type() );
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( NonlinearSolverOptionsTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
