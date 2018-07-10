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
  class AntiochBuildSpeciesTest : public CppUnit::TestCase,
                                  public GRINS::AntiochMixtureBuilderBase
  {
  public:

    CPPUNIT_TEST_SUITE( AntiochBuildSpeciesTest );

    CPPUNIT_TEST( test_ascii );

    CPPUNIT_TEST_SUITE_END();

  public:

    void test_ascii()
    {
      std::string inputfile = this->setup_ascii_input("ozone_species_data.dat");
      std::vector<std::string> species_names = this->do_work(inputfile);
      this->check_ascii(species_names);
    }

  private:

    std::vector<std::string> do_work( const std::string & inputfile )
    {
      std::stringstream inputstream;
      inputstream << inputfile;

      GetPot input(inputstream);

      std::vector<std::string> species_names;
      this->build_species_names(input,"TestMaterial",species_names);

      return species_names;
    }

    std::string setup_ascii_input( const std::string & chemical_data_filename )
    {
      std::string text = this->basic_input_setup();
      text += "chemical_data = './input_files/"+chemical_data_filename+"'\n";
      text += "species = 'O O2 O3'\n";

      return text;
    }

    std::string basic_input_setup()
    {
      std::string text = "[Materials]\n";
      text +="[./TestMaterial]\n";
      text +="[./GasMixture]\n";
      text +="[./Antioch]\n";
      return text;
    }

    void check_ascii( const std::vector<std::string> & species_names )
    {
      CPPUNIT_ASSERT_EQUAL(3,(int)species_names.size());
      CPPUNIT_ASSERT_EQUAL(std::string("O"),species_names[0]);
      CPPUNIT_ASSERT_EQUAL(std::string("O2"),species_names[1]);
      CPPUNIT_ASSERT_EQUAL(std::string("O3"),species_names[2]);
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochBuildSpeciesTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_ANTIOCH
#endif // GRINS_HAVE_CPPUNIT
