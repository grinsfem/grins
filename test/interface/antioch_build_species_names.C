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
    CPPUNIT_TEST( test_xml_air5sp );
    CPPUNIT_TEST( test_xml_air9sp );

    CPPUNIT_TEST_SUITE_END();

  public:

    void test_ascii()
    {
      std::string inputfile = this->setup_ascii_input("ozone_species_data.dat");
      std::vector<std::string> species_names = this->do_work(inputfile);
      this->check_ascii(species_names);
    }

    void test_xml_air5sp()
    {
      std::string inputfile = this->setup_xml_input("air.xml","air5sp");
      std::vector<std::string> species_names = this->do_work(inputfile);
      this->check_xml_air5sp(species_names);
    }

    void test_xml_air9sp()
    {
      std::string inputfile = this->setup_xml_input("air.xml","air9sp_CO2");
      std::vector<std::string> species_names = this->do_work(inputfile);
      this->check_xml_air9sp(species_names);
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

    std::string setup_xml_input( const std::string & chemical_data_filename, const std::string & mixture )
    {
      std::string text = this->basic_input_setup(chemical_data_filename);
      text += "gas_mixture = '"+mixture+"'\n";
      return text;
    }

    std::string setup_ascii_input( const std::string & chemical_data_filename )
    {
      std::string text = this->basic_input_setup(chemical_data_filename);
      text += "species = 'O O2 O3'\n";

      return text;
    }

    std::string basic_input_setup( const std::string & chemical_data_filename )
    {
      std::string text = "[Materials]\n";
      text +="[./TestMaterial]\n";
      text +="[./GasMixture]\n";
      text +="[./Antioch]\n";
      text += "chemical_data = './input_files/"+chemical_data_filename+"'\n";
      return text;
    }

    void check_ascii( const std::vector<std::string> & species_names )
    {
      CPPUNIT_ASSERT_EQUAL(3,(int)species_names.size());
      CPPUNIT_ASSERT_EQUAL(std::string("O"),species_names[0]);
      CPPUNIT_ASSERT_EQUAL(std::string("O2"),species_names[1]);
      CPPUNIT_ASSERT_EQUAL(std::string("O3"),species_names[2]);
    }

    void check_xml_air5sp( const std::vector<std::string> & species_names )
    {
      CPPUNIT_ASSERT_EQUAL(5,(int)species_names.size());
      CPPUNIT_ASSERT_EQUAL(std::string("N2"),species_names[0]);
      CPPUNIT_ASSERT_EQUAL(std::string("O2"),species_names[1]);
      CPPUNIT_ASSERT_EQUAL(std::string("NO"),species_names[2]);
      CPPUNIT_ASSERT_EQUAL(std::string("N"),species_names[3]);
      CPPUNIT_ASSERT_EQUAL(std::string("O"),species_names[4]);
    }

    void check_xml_air9sp( const std::vector<std::string> & species_names )
    {
      CPPUNIT_ASSERT_EQUAL(9,(int)species_names.size());
      CPPUNIT_ASSERT_EQUAL(std::string("N2"),species_names[0]);
      CPPUNIT_ASSERT_EQUAL(std::string("O2"),species_names[1]);
      CPPUNIT_ASSERT_EQUAL(std::string("NO"),species_names[2]);
      CPPUNIT_ASSERT_EQUAL(std::string("N"),species_names[3]);
      CPPUNIT_ASSERT_EQUAL(std::string("O"),species_names[4]);
      CPPUNIT_ASSERT_EQUAL(std::string("CO2"),species_names[5]);
      CPPUNIT_ASSERT_EQUAL(std::string("CO"),species_names[6]);
      CPPUNIT_ASSERT_EQUAL(std::string("CN"),species_names[7]);
      CPPUNIT_ASSERT_EQUAL(std::string("C"),species_names[8]);
    }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochBuildSpeciesTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_ANTIOCH
#endif // GRINS_HAVE_CPPUNIT
