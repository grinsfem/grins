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

// Testing
#include "grins_test_paths.h"

//GRINS
#include "grins/antioch_mixture_builder_base.h"

namespace GRINSTesting
{
  class AntiochMixtureBuilderBaseTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( AntiochMixtureBuilderBaseTest );

    CPPUNIT_TEST( test_build_chem_mix );
    CPPUNIT_TEST( test_build_reaction_set );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/antioch.in";
      _input.reset( new GetPot(input_file) );
    }

    void test_build_chem_mix()
    {
      GRINS::AntiochMixtureBuilderBase builder;

      libMesh::UniquePtr<Antioch::ChemicalMixture<libMesh::Real> > chem_mix
        = builder.build_chem_mix(*_input, "TestMaterial");

      // Just a couple of basic tests for the parsed ChemicalMixture
      CPPUNIT_ASSERT_EQUAL( 5, (int)chem_mix->n_species() );

      const std::vector<Antioch::ChemicalSpecies<libMesh::Real>*> & all_species =
        chem_mix->chemical_species();

      CPPUNIT_ASSERT_EQUAL( 5, (int)all_species.size() );
    }

    void test_build_reaction_set()
    {
      GRINS::AntiochMixtureBuilderBase builder;

      libMesh::UniquePtr<Antioch::ChemicalMixture<libMesh::Real> > chem_mix
        = builder.build_chem_mix(*_input, "TestMaterial");

      libMesh::UniquePtr<Antioch::ReactionSet<libMesh::Real> > reaction_set
        = builder.build_reaction_set(*_input, "TestMaterial", *chem_mix);

      CPPUNIT_ASSERT_EQUAL( 5, (int)reaction_set->n_species() );
      CPPUNIT_ASSERT_EQUAL( 5, (int)reaction_set->n_reactions() );
    }

  private:

    libMesh::UniquePtr<GetPot> _input;
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochMixtureBuilderBaseTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_ANTIOCH
#endif // GRINS_HAVE_CPPUNIT
