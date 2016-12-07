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
#include "system_helper.h"

// GRINS
#include "grins/builder_helper.h"

// Ignore warnings from auto_ptr in CPPUNIT_TEST_SUITE_END()
#include <libmesh/ignore_warnings.h>

namespace GRINSTesting
{
  class BuilderHelperTest : public CppUnit::TestCase,
                            public SystemHelper,
                            public GRINS::BuilderHelper // So we can test proctected methods
  {
  public:
    CPPUNIT_TEST_SUITE( BuilderHelperTest );

    CPPUNIT_TEST( test_parse_var_sections );

    CPPUNIT_TEST_SUITE_END();

  public:

    void tearDown()
    {
      this->reset_all();
    }

    void test_parse_var_sections()
    {
      std::string filename = std::string(GRINS_TEST_UNIT_INPUT_SRCDIR)+"/default_bc_builder.in";
      this->setup_multiphysics_system(filename);

      // Now we can parse the variable sections and names
      std::set<std::string> sections;
      this->parse_var_sections(*_input,sections);

      // Make sure we have the right sections
      CPPUNIT_ASSERT_EQUAL(4,(int)sections.size());
      CPPUNIT_ASSERT( sections.find("Velocity") != sections.end() );
      CPPUNIT_ASSERT( sections.find("Pressure") != sections.end() );
      CPPUNIT_ASSERT( sections.find("Temperature") != sections.end() );
      CPPUNIT_ASSERT( sections.find("SpeciesMassFractions") != sections.end() );
    }
 };

  CPPUNIT_TEST_SUITE_REGISTRATION( BuilderHelperTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
