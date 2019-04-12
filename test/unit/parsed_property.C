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

// C++
#include <string>
#include <limits>
#include <cmath>
#include <sstream>

// GRINS
#include "grins/parsed_pressure.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/point.h"

namespace GRINSTesting
{
  class ParsedPressureTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( ParsedPressureTest );

    CPPUNIT_TEST( test_parsed_pressure );

    CPPUNIT_TEST_SUITE_END();

  private:

    std::unique_ptr<GetPot> _input;

     std::string setup_input()
    {
      std::string text = "[Test]\n";
      text += "[./TestSection]\n";
      text += "pressure = 'sqrt(x^2+y^2)'";
      return text;
    }

    libMesh::Real tol()
    { return std::numeric_limits<libMesh::Real>::epsilon() * 10; }

  public:

    void test_parsed_pressure()
    {
      GRINS::ParsedPressure pressure(*_input,"Test/TestSection");
      libMesh::Point x(1.0,1.0,0.0);
      libMesh::Real value = pressure(x,1); // Not really a function of time, but let's see
      CPPUNIT_ASSERT_DOUBLES_EQUAL(std::sqrt(2.0),value,this->tol());
    }

    void setUp()
    {
      std::string input_string = this->setup_input();

      std::stringstream ss;
      ss << input_string;

      _input.reset(new GetPot(ss));
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( ParsedPressureTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
