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
  class AntiochBuilderFormatTestBase
  {
  protected:

    // Handcoding our inputfile. Common things between all formats
    std::string basic_input_setup( const std::string & chemical_data_filename )
    {
      std::string text = "[Materials]\n";
      text +="[./TestMaterial]\n";
      text +="[./GasMixture]\n";
      text +="[./Antioch]\n";
      text += "chemical_data = './input_files/"+chemical_data_filename+"'\n";
      return text;
    }

    // Specific for ASCII type parsing
    std::string setup_ascii_input( const std::string & chemical_data_filename )
    {
      std::string text = this->basic_input_setup(chemical_data_filename);
      text += "species = 'O O2 O3'\n";

      return text;
    }

    // Specific for XML type parsing
    std::string setup_xml_input( const std::string & chemical_data_filename, const std::string & mixture )
    {
      std::string text = this->basic_input_setup(chemical_data_filename);
      text += "gas_mixture = '"+mixture+"'\n";
      return text;
    }
  };

  //! Test AntiochMixtureBuilderBase::build_species_names for different parsing formats
  class AntiochBuilderFormatSpeciesNameTest : public CppUnit::TestCase,
                                              public AntiochBuilderFormatTestBase,
                                              public GRINS::AntiochMixtureBuilderBase
  {
  public:

    CPPUNIT_TEST_SUITE( AntiochBuilderFormatSpeciesNameTest );

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

  //! Test AntiochMixtureBuilderBase::build_reaction_set for different parsing formats
  class AntiochBuilderFormatKineticsTest : public CppUnit::TestCase,
                                           public AntiochBuilderFormatTestBase,
                                           public GRINS::AntiochMixtureBuilderBase
  {
  public:

    CPPUNIT_TEST_SUITE( AntiochBuilderFormatKineticsTest );

    CPPUNIT_TEST( test_ascii );
    CPPUNIT_TEST( test_xml_air5sp );
    CPPUNIT_TEST( test_xml_air9sp );

    CPPUNIT_TEST_SUITE_END();

  public:

    void test_ascii()
    {
      std::string inputfile = this->setup_ascii_kinetics_input("ozone_species_data.dat",
                                                               "ozone.xml");
      std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > mixture;
      std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > reaction_set;

      this->do_work( inputfile, mixture, reaction_set );

      this->check_ozone_kinetics(*mixture,*reaction_set);
    }

    void test_xml_air5sp()
    {
      std::string inputfile = this->setup_xml_input("air.xml","air5sp");
      std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > mixture;
      std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > reaction_set;

      this->do_work( inputfile, mixture, reaction_set );

      this->check_air5sp_kinetics(*mixture,*reaction_set);
    }

    void test_xml_air9sp()
    {
      std::string inputfile = this->setup_xml_input("air.xml","air9sp_CO2");
      std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > mixture;
      std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > reaction_set;

      this->do_work( inputfile, mixture, reaction_set );

      this->check_air9sp_kinetics(*mixture,*reaction_set);
    }

  private:

    void do_work( const std::string & inputfile,
                  std::unique_ptr<Antioch::ChemicalMixture<libMesh::Real> > & mixture,
                  std::unique_ptr<Antioch::ReactionSet<libMesh::Real> > & reaction_set )
    {
      std::stringstream inputstream;
      inputstream << inputfile;

      GetPot input(inputstream);

      mixture = this->build_chem_mix(input,"TestMaterial");
      reaction_set = this->build_reaction_set(input,"TestMaterial",*mixture);
    }

    std::string setup_ascii_kinetics_input( const std::string & chemical_data_filename,
                                            const std::string & kinetics_filename )
    {
      std::string text = this->setup_ascii_input(chemical_data_filename);
      text += "kinetics_data = './input_files/"+kinetics_filename+"'\n";

      return text;
    }

    void check_ozone_kinetics( const Antioch::ChemicalMixture<libMesh::Real> & mixture,
                               const Antioch::ReactionSet<libMesh::Real> & reaction_set )
    {
      CPPUNIT_ASSERT_EQUAL(3,(int)reaction_set.n_reactions());

      // Check first reaction
      {
        const Antioch::Reaction<libMesh::Real> & modarr_threebody = reaction_set.reaction(0);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, modarr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, modarr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(modarr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            libMesh::Real efficiency = modarr_threebody.get_efficiency(s);

            if( species_name == std::string("O3") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.92, efficiency, this->tol());
            else if( species_name == std::string("O2") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.94, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.13, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }


      // Check second reaction
      {
        const Antioch::Reaction<libMesh::Real> & arr_threebody = reaction_set.reaction(1);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, arr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::ARRHENIUS, arr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(arr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            libMesh::Real efficiency = arr_threebody.get_efficiency(s);

            if( species_name == std::string("O3") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.92, efficiency, this->tol());
            else if( species_name == std::string("O2") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(0.94, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.13, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }

      // Check third reaction
      {
        const Antioch::Reaction<libMesh::Real> & arr_rev = reaction_set.reaction(2);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,arr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::ARRHENIUS, arr_rev.kinetics_model() );
        CPPUNIT_ASSERT(arr_rev.reversible());
      }
    }


    void check_air5sp_kinetics( const Antioch::ChemicalMixture<libMesh::Real> & mixture,
                                const Antioch::ReactionSet<libMesh::Real> & reaction_set )
    {
      // This file actually has a lot more reactions, so we're additionally
      // testing that we're stripping off reactions that include species
      // that are not in the mixture that formed the reaction set.
      CPPUNIT_ASSERT_EQUAL(5,(int)reaction_set.n_reactions());

      this->check_air5sp_reactions_only(mixture,reaction_set);
    }

    void check_air9sp_kinetics( const Antioch::ChemicalMixture<libMesh::Real> & mixture,
                                const Antioch::ReactionSet<libMesh::Real> & reaction_set )
    {
      // This file actually has a lot more reactions, so we're additionally
      // testing that we're stripping off reactions that include species
      // that are not in the mixture that formed the reaction set.
      CPPUNIT_ASSERT_EQUAL(15,(int)reaction_set.n_reactions());

      this->check_air9sp_reactions_only(mixture,reaction_set);
    }

    void check_air9sp_reactions_only( const Antioch::ChemicalMixture<libMesh::Real> & mixture,
                                      const Antioch::ReactionSet<libMesh::Real> & reaction_set )
    {
      // These reactions should be the same as the air5sp case
      this->check_air5sp_reactions_only(mixture,reaction_set);

      // Now onto the other 10 reactions...

      // Reaction 7 in the air.xml file, reaction 5 in our ReactionSet
      {
        const Antioch::Reaction<libMesh::Real> & arr_threebody = reaction_set.reaction(5);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, arr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::ARRHENIUS, arr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(arr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            libMesh::Real efficiency = arr_threebody.get_efficiency(s);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }

      // Reaction 10 in the air.xml file, reaction 6 in our ReactionSet
      {
        const Antioch::Reaction<libMesh::Real> & arr_rev = reaction_set.reaction(6);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,arr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, arr_rev.kinetics_model() );
        CPPUNIT_ASSERT(arr_rev.reversible());
      }

      // Reaction 11 in the air.xml file, reaction 7 in our ReactionSet
      {
        const Antioch::Reaction<libMesh::Real> & arr_rev = reaction_set.reaction(7);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,arr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::ARRHENIUS, arr_rev.kinetics_model() );
        CPPUNIT_ASSERT(arr_rev.reversible());
      }

      // Reaction 12 in the air.xml file, reaction 8 in our ReactionSet
      {
        const Antioch::Reaction<libMesh::Real> & arr_rev = reaction_set.reaction(8);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,arr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, arr_rev.kinetics_model() );
        CPPUNIT_ASSERT(arr_rev.reversible());
      }

      // Reaction 13 in the air.xml file, reaction 9 in our ReactionSet
      {
        const Antioch::Reaction<libMesh::Real> & arr_rev = reaction_set.reaction(9);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,arr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, arr_rev.kinetics_model() );
        CPPUNIT_ASSERT(arr_rev.reversible());
      }

      // Reaction 26 in the air.xml file, reaction 10 in our ReactionSet
      {
        const Antioch::Reaction<libMesh::Real> & arr_rev = reaction_set.reaction(10);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,arr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, arr_rev.kinetics_model() );
        CPPUNIT_ASSERT(arr_rev.reversible());
      }

      // Reaction 27 in the air.xml file, reaction 11 in our ReactionSet
      {
        const Antioch::Reaction<libMesh::Real> & arr_rev = reaction_set.reaction(11);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,arr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, arr_rev.kinetics_model() );
        CPPUNIT_ASSERT(arr_rev.reversible());
      }

      // Reaction 28 in the air.xml file, reaction 12 in our ReactionSet
      {
        const Antioch::Reaction<libMesh::Real> & arr_rev = reaction_set.reaction(12);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,arr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::ARRHENIUS, arr_rev.kinetics_model() );
        CPPUNIT_ASSERT(arr_rev.reversible());
      }

      // Reaction 29 in the air.xml file, reaction 13 in our ReactionSet
      {
        const Antioch::Reaction<libMesh::Real> & arr_threebody = reaction_set.reaction(13);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, arr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, arr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(arr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            libMesh::Real efficiency = arr_threebody.get_efficiency(s);

            if( species_name == std::string("N") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(2.029, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(2.029, efficiency, this->tol());
            else if( species_name == std::string("C") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(2.029, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }

      // Reaction 30 in the air.xml file, reaction 14 in our ReactionSet
      {
        const Antioch::Reaction<libMesh::Real> & arr_threebody = reaction_set.reaction(14);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, arr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, arr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(arr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            libMesh::Real efficiency = arr_threebody.get_efficiency(s);

            if( species_name == std::string("N") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.478, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.478, efficiency, this->tol());
            else if( species_name == std::string("C") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.478, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }

    }


    void check_air5sp_reactions_only( const Antioch::ChemicalMixture<libMesh::Real> & mixture,
                                      const Antioch::ReactionSet<libMesh::Real> & reaction_set )
    {
      // Check first reaction
      {
        const Antioch::Reaction<libMesh::Real> & modarr_threebody = reaction_set.reaction(0);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, modarr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, modarr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(modarr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            libMesh::Real efficiency = modarr_threebody.get_efficiency(s);

            if( species_name == std::string("N") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(4.2857, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(4.2857, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }

      // Check second reaction
      {
        const Antioch::Reaction<libMesh::Real> & modarr_threebody = reaction_set.reaction(1);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, modarr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, modarr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(modarr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            libMesh::Real efficiency = modarr_threebody.get_efficiency(s);

            if( species_name == std::string("N") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }

      // Check third reaction
      {
        const Antioch::Reaction<libMesh::Real> & arr_threebody = reaction_set.reaction(2);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::THREE_BODY, arr_threebody.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::ARRHENIUS, arr_threebody.kinetics_model() );
        CPPUNIT_ASSERT(arr_threebody.reversible());

        // Check threebody efficiencies
        for( unsigned int s = 0; s < reaction_set.n_species(); s++ )
          {
            std::string species_name = mixture.chemical_species()[s]->species();
            libMesh::Real efficiency = arr_threebody.get_efficiency(s);

            if( species_name == std::string("N") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(22.0, efficiency, this->tol());
            else if( species_name == std::string("O") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(22.0, efficiency, this->tol());
            else if( species_name == std::string("NO") )
              CPPUNIT_ASSERT_DOUBLES_EQUAL(22.0, efficiency, this->tol());
            else
              // Default should be 1.0
              CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, efficiency, this->tol());
          }
      }


      // Check fourth reaction
      {
        const Antioch::Reaction<libMesh::Real> & modarr_rev = reaction_set.reaction(3);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,modarr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::KOOIJ, modarr_rev.kinetics_model() );
        CPPUNIT_ASSERT(modarr_rev.reversible());
      }

      // Check fifth reaction
      {
        const Antioch::Reaction<libMesh::Real> & arr_rev = reaction_set.reaction(4);
        CPPUNIT_ASSERT_EQUAL(Antioch::ReactionType::ELEMENTARY ,arr_rev.type() );
        CPPUNIT_ASSERT_EQUAL(Antioch::KineticsModel::ARRHENIUS, arr_rev.kinetics_model() );
        CPPUNIT_ASSERT(arr_rev.reversible());
      }
    }

    libMesh::Real tol()
    { return std::numeric_limits<libMesh::Real>::epsilon() * 10; }

  };

  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochBuilderFormatSpeciesNameTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochBuilderFormatKineticsTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_ANTIOCH
#endif // GRINS_HAVE_CPPUNIT
