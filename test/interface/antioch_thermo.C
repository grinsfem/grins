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
#ifdef GRINS_HAVE_ANTIOCH

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

#include "grins_test_paths.h"
#include "air_nasa_poly_base.h"
#include "antioch_test_base.h"
#include "thermochem_test_common.h"
#include "testing_utils.h"

// GRINS
#include "grins/antioch_evaluator.h"

// C++
#include <sstream>

namespace GRINSTesting
{
  class AntiochCEAThermoTest : public CppUnit::TestCase,
                               public AntiochTestBase,
                               public AirNASA9TestBase
  {
  public:
    CPPUNIT_TEST_SUITE( AntiochCEAThermoTest );

    CPPUNIT_TEST( test_cp );
    CPPUNIT_TEST( test_hs );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/antioch.in";
      this->init_antioch(input_file, "TestMaterial");

      this->check_indices();
    }

    void check_indices()
    {
      CPPUNIT_ASSERT_EQUAL( _N2_idx, _antioch_mixture->species_index("N2") );
      CPPUNIT_ASSERT_EQUAL( _O2_idx,  _antioch_mixture->species_index("O2") );
      CPPUNIT_ASSERT_EQUAL( _N_idx,  _antioch_mixture->species_index("N") );
      CPPUNIT_ASSERT_EQUAL( _O_idx,  _antioch_mixture->species_index("O") );
      CPPUNIT_ASSERT_EQUAL( _NO_idx, _antioch_mixture->species_index("NO") );
    }

    void test_cp()
    {
      this->test_cp_common<GRINS::AntiochMixture,GRINS::AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> > >
        ( *_antioch_mixture, (*this), TestingUtils::epsilon()*100 );
    }

    void test_hs()
    {
      this->test_h_common<GRINS::AntiochMixture,GRINS::AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> > >
        ( *_antioch_mixture, (*this), TestingUtils::epsilon()*100 );
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochCEAThermoTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_ANTIOCH
#endif // GRINS_HAVE_CPPUNIT
