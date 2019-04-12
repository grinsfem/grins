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

      libMesh::Real tolerance = 1.0e-9;

      // test getting arbitrary data values
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0,hitran.isotopologue(0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3,hitran.isotopologue(20),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3682.70083,hitran.nu0(0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3682.761088,hitran.nu0(20),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.062e-30,hitran.sw(0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.48e-30,hitran.sw(20),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0685,hitran.gamma_air(0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0679,hitran.gamma_air(20),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.084,hitran.gamma_self(0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08,hitran.gamma_self(20),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3253.949,hitran.elower(0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1349.8942,hitran.elower(20),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.76,hitran.n_air(0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(0.75,hitran.n_air(20),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.008269,hitran.delta_air(0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.004588,hitran.delta_air(20),tolerance);


      // these T values are explicitly given in the partition sum data
      CPPUNIT_ASSERT_DOUBLES_EQUAL(279.609573308,hitran.partition_function(290.02,0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(563.425588693,hitran.partition_function(290.02,1),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(593.852717624,hitran.partition_function(290.02,2),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3461.78001223,hitran.partition_function(290.02,3),tolerance);

      // check T value on the ends
      CPPUNIT_ASSERT_DOUBLES_EQUAL(279.585269287,hitran.partition_function(290.0,0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(563.375898926,hitran.partition_function(290.0,1),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(593.800881348,hitran.partition_function(290.0,2),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3461.47846875,hitran.partition_function(290.0,3),tolerance);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(304.559997559,hitran.partition_function(310.0,0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(614.489990234,hitran.partition_function(310.0,1),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(647.090026855,hitran.partition_function(310.0,2),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3771.39990234,hitran.partition_function(310.0,3),tolerance);


      // partition function values at the reference temperature T=296K
      CPPUNIT_ASSERT_DOUBLES_EQUAL(286.935573058,hitran.partition_function(296,0),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(578.408361459,hitran.partition_function(296,1),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(609.479752969,hitran.partition_function(296,2),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(3552.67876127,hitran.partition_function(296,3),tolerance);

      // this T value is not is the partition sum data, and must be (linearly) interpolated
      CPPUNIT_ASSERT_DOUBLES_EQUAL(593.813840240,hitran.partition_function(290.005,2),tolerance);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(291.9025545935294,hitran.partition_function(300+1.0e-6,0),1.0e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(291.9137901096000,hitran.partition_function(300.009,0),1.0e-13);

    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( HITRANtest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
