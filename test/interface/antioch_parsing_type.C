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
#ifdef GRINS_HAVE_ANTIOCH

#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>

//GRINS
#include "grins/antioch_mixture_builder_base.h"

namespace GRINSTesting
{
  class AntiochParsingTypeTest : public CppUnit::TestCase,
                                 public GRINS::AntiochMixtureBuilderBase
  {
  public:

    CPPUNIT_TEST_SUITE( AntiochParsingTypeTest );

    CPPUNIT_TEST( test_ascii );
    CPPUNIT_TEST( test_xml );
    CPPUNIT_TEST( test_chemkin );

    CPPUNIT_TEST_SUITE_END();

  public:

    void test_ascii()
    {
      std::string suffix(".dat");
      this->run_test(suffix,Antioch::ASCII);
    }

    void test_xml()
    {
      std::string suffix(".xml");
      this->run_test(suffix,Antioch::XML);
    }

    void test_chemkin()
    {
      std::string suffix(".chemkin");
      this->run_test(suffix,Antioch::CHEMKIN);
    }

  private:

    void run_test( const std::string & suffix, Antioch::ParsingType exact_parsing_type )
    {
      GetPot input;
      this->setup_inputfile(suffix,input);

      Antioch::ParsingType parsing_type =
        this->get_antioch_parsing_type(input,std::string("TestMaterial"));

      CPPUNIT_ASSERT_EQUAL(exact_parsing_type,parsing_type);
    }

    void setup_inputfile( const std::string & suffix, GetPot & input )
    {
      std::string filename = "testfile"+suffix;

      std::string text = "[Materials]\n";
      text +="[./TestMaterial]\n";
      text +="[./GasMixture]\n";
      text +="[./Antioch]\n";
      text += "chemical_data = '"+filename+"'\n";

      std::stringstream inputfile;
      inputfile << text;
      input.parse_input_stream(inputfile);
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochParsingTypeTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_ANTIOCH
#endif // GRINS_HAVE_CPPUNIT
