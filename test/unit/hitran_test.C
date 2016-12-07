//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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

// Testing headers
#include "test_comm.h"
#include "grins_test_paths.h"

// GRINS
#include "grins/hitran.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class HITRANtest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( HITRANtest );

    CPPUNIT_TEST( parse_from_file );

    CPPUNIT_TEST_SUITE_END();

  public:

    void parse_from_file()
    {
      std::string data_file = std::string(GRINS_TEST_SRCDIR)+"/test_data/CO2_data.dat";
      std::string partition_file = std::string(GRINS_TEST_SRCDIR)+"/test_data/CO2_partition_function.dat";
      
      libMesh::Real T_min = 290;
      libMesh::Real T_max = 310;
      libMesh::Real T_step = 0.01;

      GRINS::HITRAN hitran(data_file,partition_file,T_min,T_max,T_step);

      // test getting arbitrary data values
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0,hitran.isotopologue(0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3,hitran.isotopologue(20),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3682.70083,hitran.nu0(0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3682.761088,hitran.nu0(20),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.062e-30,hitran.sw(0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.48e-30,hitran.sw(20),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0685,hitran.gamma_air(0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0679,hitran.gamma_air(20),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.084,hitran.gamma_self(0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08,hitran.gamma_self(20),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3253.949,hitran.elower(0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1349.8942,hitran.elower(20),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.76,hitran.n_air(0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.75,hitran.n_air(20),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.008269,hitran.delta_air(0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.004588,hitran.delta_air(20),libMesh::TOLERANCE);
      

      // these T values are explicitly given in the partition sum data
      CPPUNIT_ASSERT_DOUBLES_EQUAL(279.609573308,hitran.partition_function(290.02,0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(563.425588693,hitran.partition_function(290.02,1),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(593.852717624,hitran.partition_function(290.02,2),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3461.78001223,hitran.partition_function(290.02,3),libMesh::TOLERANCE);
      
      // partition function values at the reference temperature T=296K
      CPPUNIT_ASSERT_DOUBLES_EQUAL(286.93557306,hitran.partition_function(296,0),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(578.40836146,hitran.partition_function(296,1),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(609.47975297,hitran.partition_function(296,2),libMesh::TOLERANCE);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3552.67876127,hitran.partition_function(296,3),libMesh::TOLERANCE);

      // this T value is not is the partition sum data, and must be (linearly) interpolated
      CPPUNIT_ASSERT_DOUBLES_EQUAL(593.81384024,hitran.partition_function(290.005,2),libMesh::TOLERANCE);

    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( HITRANtest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
