//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

// libMesh
#include "libmesh/parsed_function.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class LaserAbsorptionTest : public CppUnit::TestCase,
                              public SpectroscopicTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( LaserAbsorptionTest );
    CPPUNIT_TEST( two_dimension_mesh );
    CPPUNIT_TEST( three_dimension_mesh );
    CPPUNIT_TEST_SUITE_END();

  public:

    void tearDown()
    {
      // Clear out the VariableWarehouse so it doesn't interfere with other tests.
      GRINS::GRINSPrivate::VariableWarehouse::clear();
    }

    //! 10x10 mesh of QUAD4 with simple, non-moving flow. Hence, orientation and position of the optical paths
    //! should not change the QoI value
    void two_dimension_mesh()
    {
      libMesh::Real calc_answer = 4.7959591239712063e-01;

      // left to right
      std::stringstream ss;
      ss << this->laser_string_2D("laser_absorption","LaserAbsorption","0.0 0.0501","0.0 0.0500","0.0 0.0499",0.0,10,10);
      this->run_test(ss,calc_answer);

      this->tearDown();

      // bottom to top
      std::stringstream ss2;
      ss2 << this->laser_string_2D("laser_absorption","LaserAbsorption","0.0501 0.0","0.0500 0.0","0.0499 0.0",1.57079632679,10,10);
      this->run_test(ss2,calc_answer);
    }

    //! 3x3x3 mesh of HEX8 with simple, non-moving flow. Hence, orientation and position of the optical paths
    //! should not change the QoI value
    void three_dimension_mesh()
    {
      libMesh::Real calc_answer = 4.7959591239712063e-01;

      // left to right
      std::stringstream ss;
      ss << this->laser_string_3D("laser_absorption","LaserAbsorption","0.0500 0.0501 0.0","0.0500 0.0500 0.0","0.0501 0.0500 0.0",0.0,0.0,3,3,3);
      this->run_test(ss,calc_answer);

      this->tearDown();

      // bottom to top
      std::stringstream ss2;
      ss2 << this->laser_string_3D("laser_absorption","LaserAbsorption","0.0 0.0500 0.0501","0.0 0.0500 0.0500","0.0 0.0501 0.0500",0.0,1.57079632679,3,3,3);
      this->run_test(ss2,calc_answer);

      this->tearDown();

      // front to back
      std::stringstream ss3;
      ss3 << this->laser_string_3D("laser_absorption","LaserAbsorption","0.0500 0.0 0.0501","0.0500 0.0 0.0500","0.0501 0.0 0.0500",1.57079632679,1.57079632679,3,3,3);
      this->run_test(ss3,calc_answer);
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( LaserAbsorptionTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_ANTIOCH
#endif // GRINS_HAVE_CPPUNIT
