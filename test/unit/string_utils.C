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

#include <string>
#include <vector>
#include <iostream>
#include <limits>

#include "grins/string_utils.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class StringUtilitiesTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( StringUtilitiesTest );

    CPPUNIT_TEST( test_split_string );
    CPPUNIT_TEST( test_string_to_T );
    CPPUNIT_TEST( test_T_to_string );

    CPPUNIT_TEST_SUITE_END();

  public:

    void test_split_string()
    {
      {
        std::string str_1("N->N2");
        std::vector<std::string> test_1_split_exact(2);
        test_1_split_exact[0] = std::string("N");
        test_1_split_exact[1] = std::string("N2");

        std::vector<std::string> str_1_split;
        GRINS::StringUtilities::split_string( str_1, "->", str_1_split);
        this->test_string( str_1_split, test_1_split_exact );
      }

      {
        std::string str_2("N+C(s)->CN");
        std::vector<std::string> test_2_split_exact(2);
        test_2_split_exact[0] = std::string("N+C(s)");
        test_2_split_exact[1] = std::string("CN");

        std::vector<std::string> str_2_split;
        GRINS::StringUtilities::split_string( str_2, "->", str_2_split);
        this->test_string( str_2_split, test_2_split_exact );
      }

      {
        std::string str_3("u:v:w:T:p:w_N:w_N2:p0");
        std::vector<std::string> test_3_split_exact(8);
        test_3_split_exact[0] = std::string("u");
        test_3_split_exact[1] = std::string("v");
        test_3_split_exact[2] = std::string("w");
        test_3_split_exact[3] = std::string("T");
        test_3_split_exact[4] = std::string("p");
        test_3_split_exact[5] = std::string("w_N");
        test_3_split_exact[6] = std::string("w_N2");
        test_3_split_exact[7] = std::string("p0");

        std::vector<std::string> str_3_split;
        GRINS::StringUtilities::split_string( str_3, ":", str_3_split);
        this->test_string( str_3_split, test_3_split_exact );
      }

      {
        std::string str_4("u v w T p w_N w_N2 p0");
        std::vector<std::string> test_4_split_exact(8);
        test_4_split_exact[0] = std::string("u");
        test_4_split_exact[1] = std::string("v");
        test_4_split_exact[2] = std::string("w");
        test_4_split_exact[3] = std::string("T");
        test_4_split_exact[4] = std::string("p");
        test_4_split_exact[5] = std::string("w_N");
        test_4_split_exact[6] = std::string("w_N2");
        test_4_split_exact[7] = std::string("p0");

        std::vector<std::string> str_4_split;
        GRINS::StringUtilities::split_string( str_4, " ", str_4_split);
        this->test_string( str_4_split, test_4_split_exact );
      }
    }

    void test_string_to_T()
    {
      std::string one = "1";
      int ione = GRINS::StringUtilities::string_to_T<int>(one);
      unsigned int uione = GRINS::StringUtilities::string_to_T<unsigned int>(one);
      CPPUNIT_ASSERT_EQUAL(1,ione);
      CPPUNIT_ASSERT_EQUAL((unsigned int)1,uione);

      std::string tenp1 = "10.1";
      double dtenp1 = GRINS::StringUtilities::string_to_T<double>(tenp1);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(10.1,
                                   dtenp1,
                                   std::numeric_limits<double>::epsilon());
    }

    void test_T_to_string()
    {
      std::string one_exact = "1";
      std::string tenp1_exact = "10.1";

      {
        std::string one_test =  GRINS::StringUtilities::T_to_string<int>(1);
        CPPUNIT_ASSERT_EQUAL(one_exact,one_test);
      }

      {
        std::string one_test =  GRINS::StringUtilities::T_to_string<unsigned int>(1);
        CPPUNIT_ASSERT_EQUAL(one_exact,one_test);
      }

      {
        std::string tenp1_test =  GRINS::StringUtilities::T_to_string<double>(10.1);
        CPPUNIT_ASSERT_EQUAL(tenp1_exact,tenp1_test);
      }
    }

  private:

    void test_string( const std::vector<std::string>& test,
                      const std::vector<std::string>& exact )
    {
      CPPUNIT_ASSERT_EQUAL(test.size(), exact.size() );

      for( unsigned int s = 0; s < test.size(); s++ )
        CPPUNIT_ASSERT_EQUAL( test[s], exact[s] );
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( StringUtilitiesTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
