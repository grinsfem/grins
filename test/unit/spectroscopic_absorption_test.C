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

#ifdef GRINS_HAVE_ANTIOCH

#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include "test_comm.h"
#include "grins_test_paths.h"

#include "spectroscopic_test_base.h"

// GRINS
#include "grins/math_constants.h"
#include "grins/grins_enums.h"
#include "grins/simulation_builder.h"
#include "grins/simulation.h"
#include "grins/variable_warehouse.h"
#include "grins/single_variable.h"
#include "grins/multicomponent_variable.h"
#include "grins/chemistry_builder.h"
#include "grins/antioch_chemistry.h"

// libMesh
#include "libmesh/parsed_function.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class SpectroscopicAbsorptionTest : public CppUnit::TestCase,
                                      public SpectroscopicTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( SpectroscopicAbsorptionTest );
    CPPUNIT_TEST( single_elem_mesh );
    CPPUNIT_TEST( multi_elem_mesh );
    CPPUNIT_TEST( param_derivs );
    CPPUNIT_TEST( elem_qoi_derivatives );
    CPPUNIT_TEST( test_assemble_qoi_derivatives );
    CPPUNIT_TEST_SUITE_END();

  public:

    void tearDown()
    {
      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

    //! Single QUAD4 elem, uniform T,P,Y
    void single_elem_mesh()
    {
      std::stringstream ss;
      ss << this->spectroscopic_string("spectroscopic_absorption","SpectroscopicAbsorption",1,1);

      libMesh::Real calc_answer = 1.0 - 0.520403290868787;

      this->run_test(ss,calc_answer);
    }

    //! 10x10 mesh, uniform T,P,Y
    void multi_elem_mesh()
    {
      std::stringstream ss;
      ss << this->spectroscopic_string("spectroscopic_absorption","SpectroscopicAbsorption",3,3);

      libMesh::Real calc_answer = 1.0 - 0.520403290868787; // same physical conditions as single_elem_mesh, so answer should not change

      this->run_test(ss,calc_answer);
    }

    void param_derivs()
    {
      std::stringstream ss;
      ss << this->spectroscopic_string("spectroscopic_absorption","SpectroscopicAbsorption",1,1);

      this->param_deriv_test(ss);
    }

    void elem_qoi_derivatives()
    {
      std::stringstream ss;
      ss << this->spectroscopic_string("spectroscopic_absorption","SpectroscopicAbsorption",1,1);

      this->elem_qoi_derivative_test(ss);
    }

    void test_assemble_qoi_derivatives()
    {
      std::stringstream ss;
      ss << this->spectroscopic_string("spectroscopic_absorption","SpectroscopicAbsorption",1,1);

      this->assemble_qoi_derivative_test(ss);
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( SpectroscopicAbsorptionTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_ANTIOCH
#endif // GRINS_HAVE_CPPUNIT
