//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

// Testing headers
#include "test_comm.h"
#include "grins_test_paths.h"
#include "system_helper.h"

// GRINS
#include "grins/default_bc_builder.h"

namespace GRINSTesting
{
  class DefaultBCBuilderTest : public CppUnit::TestCase,
                               public GRINS::DefaultBCBuilder // So we can test proctected methods
  {
  public:
    CPPUNIT_TEST_SUITE( DefaultBCBuilderTest );

    CPPUNIT_TEST_SUITE_END();

  public:

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( DefaultBCBuilderTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
