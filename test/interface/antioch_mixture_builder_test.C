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
#include "grins/antioch_constant_transport_mixture_builder.h"
#include "grins/antioch_mixture_averaged_transport_mixture_builder.h"

namespace GRINSTesting
{
  class AntiochMixtureBuilderBaseTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( AntiochMixtureBuilderBaseTest );

    CPPUNIT_TEST( test_build_chem_mix );
    CPPUNIT_TEST( test_build_reaction_set );
    CPPUNIT_TEST( test_build_nasa_thermo_mix_cea );
    CPPUNIT_TEST( test_parse_min_T );
    CPPUNIT_TEST( test_parse_clip_negative_rho );

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

    void test_build_nasa_thermo_mix_cea()
    {
      GRINS::AntiochMixtureBuilderBase builder;

      libMesh::UniquePtr<Antioch::ChemicalMixture<libMesh::Real> > chem_mix
        = builder.build_chem_mix(*_input, "TestMaterial");

      libMesh::UniquePtr<Antioch::NASAThermoMixture<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> > >
        nasa_mix =
        builder.build_nasa_thermo_mix<Antioch::CEACurveFit<libMesh::Real> >(*_input, "TestMaterial", *chem_mix);

      CPPUNIT_ASSERT(nasa_mix->check());
    }

    void test_parse_min_T()
    {
      GRINS::AntiochMixtureBuilderBase builder;

      libMesh::Real min_T_exact = 114.15;
      libMesh::Real min_T = builder.parse_min_T( *_input, "TestMaterial" );

      CPPUNIT_ASSERT_DOUBLES_EQUAL( min_T_exact, min_T, std::numeric_limits<libMesh::Real>::epsilon() );
    }

    void test_parse_clip_negative_rho()
    {
      GRINS::AntiochMixtureBuilderBase builder;

      bool clip_negative_rho_exact = true;
      bool clip_negative_rho = builder.parse_clip_negative_rho( *_input, "TestMaterial" );

      CPPUNIT_ASSERT_EQUAL( clip_negative_rho_exact, clip_negative_rho );
    }

  private:

    libMesh::UniquePtr<GetPot> _input;
  };


  class AntiochConstantTransportMixtureBuilderTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( AntiochConstantTransportMixtureBuilderTest );

    CPPUNIT_TEST( test_build_constant_viscosity );
    CPPUNIT_TEST( test_build_constant_conductivity );
    CPPUNIT_TEST( test_build_constant_prandtl_conductivity );
    CPPUNIT_TEST( test_build_constant_lewis_diff );
    CPPUNIT_TEST( test_build_cea_constcond_mix );
    CPPUNIT_TEST( test_build_cea_constpr_mix );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/antioch.in";
      _input.reset( new GetPot(input_file) );
    }

    void test_build_constant_viscosity()
    {
      GRINS::AntiochConstantTransportMixtureBuilder builder;

      libMesh::UniquePtr<GRINS::ConstantViscosity> visc =
        builder.build_constant_viscosity( *_input, "TestMaterial" );

      libMesh::Real mu = (*visc)();

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.14, mu, std::numeric_limits<libMesh::Real>::epsilon() );
    }

    void test_build_constant_conductivity()
    {
      GRINS::AntiochConstantTransportMixtureBuilder builder;

      libMesh::UniquePtr<GRINS::ConstantConductivity> conductivity =
        builder.build_constant_conductivity<GRINS::ConstantConductivity>( *_input, "TestMaterial" );

      libMesh::Real k = (*conductivity)();

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.71, k, std::numeric_limits<libMesh::Real>::epsilon() );
    }

    void test_build_constant_prandtl_conductivity()
    {
      GRINS::AntiochConstantTransportMixtureBuilder builder;

      libMesh::UniquePtr<GRINS::ConstantPrandtlConductivity> conductivity =
        builder.build_constant_conductivity<GRINS::ConstantPrandtlConductivity>( *_input, "TestMaterial" );

      libMesh::Real mu = 1.2;
      libMesh::Real cp = 2.3;
      libMesh::Real Pr = 0.7;
      libMesh::Real k = (*conductivity)( mu, cp );
      libMesh::Real k_exact = mu*cp/Pr;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( k_exact, k, std::numeric_limits<libMesh::Real>::epsilon() );
    }

    void test_build_constant_lewis_diff()
    {
      GRINS::AntiochConstantTransportMixtureBuilder builder;

      libMesh::UniquePtr<Antioch::ConstantLewisDiffusivity<libMesh::Real> > diff =
        builder.build_constant_lewis_diff( *_input, "TestMaterial" );

      libMesh::Real rho = 1.2;
      libMesh::Real cp = 2.3;
      libMesh::Real k = 3.4;
      libMesh::Real Le = 1.4;

      libMesh::Real D = diff->D(rho,cp,k);

      libMesh::Real D_exact = k/rho/cp/Le;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( D_exact, D, std::numeric_limits<libMesh::Real>::epsilon() );
    }

    void test_build_cea_constcond_mix()
    {
      GRINS::AntiochConstantTransportMixtureBuilder builder;

      libMesh::UniquePtr<GRINS::AntiochConstantTransportMixture<Antioch::CEACurveFit<libMesh::Real>,
                                                                GRINS::ConstantConductivity> >
        mixture = builder.build_mixture<Antioch::CEACurveFit<libMesh::Real>,GRINS::ConstantConductivity>
        (*_input, "TestMaterial");

      libMesh::Real mu = mixture->mu();
      libMesh::Real k = (mixture->conductivity())();

      libMesh::Real rho = 1.2;
      libMesh::Real cp = 2.3;

      libMesh::Real Le = 1.4;

      libMesh::Real D = mixture->diffusivity().D(rho,cp,k);

      libMesh::Real D_exact = k/rho/cp/Le;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.14, mu, std::numeric_limits<libMesh::Real>::epsilon() );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( 2.71, k, std::numeric_limits<libMesh::Real>::epsilon() );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( D_exact, D, std::numeric_limits<libMesh::Real>::epsilon() );
    }

    void test_build_cea_constpr_mix()
    {
      GRINS::AntiochConstantTransportMixtureBuilder builder;

      libMesh::UniquePtr<GRINS::AntiochConstantTransportMixture<Antioch::CEACurveFit<libMesh::Real>,
                                                                GRINS::ConstantPrandtlConductivity> >
        mixture = builder.build_mixture<Antioch::CEACurveFit<libMesh::Real>,GRINS::ConstantPrandtlConductivity>
        (*_input, "TestMaterial");

      libMesh::Real mu = mixture->mu();
      libMesh::Real cp = 2.3;
      libMesh::Real Pr = 0.7;
      libMesh::Real k = (mixture->conductivity())(mu, cp);
      libMesh::Real k_exact = mu*cp/Pr;

      libMesh::Real rho = 1.2;

      libMesh::Real Le = 1.4;

      libMesh::Real D = mixture->diffusivity().D(rho,cp,k);

      libMesh::Real D_exact = k/rho/cp/Le;

      CPPUNIT_ASSERT_DOUBLES_EQUAL( 3.14, mu, std::numeric_limits<libMesh::Real>::epsilon() );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( k_exact, k, std::numeric_limits<libMesh::Real>::epsilon() );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( D_exact, D, std::numeric_limits<libMesh::Real>::epsilon() );
    }

  private:

    libMesh::UniquePtr<GetPot> _input;
  };




  class AntiochMixtureAveragedTransportMixtureBuilderTest : public CppUnit::TestCase
  {
  public:
    CPPUNIT_TEST_SUITE( AntiochMixtureAveragedTransportMixtureBuilderTest );

    CPPUNIT_TEST( test_build_cea_sutherland_eucken_constlewis_mix );
    CPPUNIT_TEST( test_build_cea_statmech_sutherland_eucken_constlewis_mix );

#ifdef ANTIOCH_HAVE_GSL
    CPPUNIT_TEST( test_build_cea_kinetic_theory_mix );
#endif // ANTIOCH_HAVE_GSL

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/antioch.in";
      _input.reset( new GetPot(input_file) );
    }

    void test_build_cea_sutherland_eucken_constlewis_mix()
    {
      GRINS::AntiochMixtureAveragedTransportMixtureBuilder builder;

      libMesh::UniquePtr<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::CEACurveFit<libMesh::Real>,
                                                                       Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real>,
                                                                       Antioch::SutherlandViscosity<libMesh::Real>,
                                                                       Antioch::EuckenThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real> > ,
                                                                       Antioch::ConstantLewisDiffusivity<libMesh::Real> > >
        mixture = builder.build_mixture<Antioch::CEACurveFit<libMesh::Real>,
                                        Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real>,
                                        Antioch::SutherlandViscosity<libMesh::Real>,
                                        Antioch::EuckenThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real> >,
                                        Antioch::ConstantLewisDiffusivity<libMesh::Real> >( *_input, "TestMaterial" );

      CPPUNIT_ASSERT(mixture);
    }

    void test_build_cea_statmech_sutherland_eucken_constlewis_mix()
    {
      GRINS::AntiochMixtureAveragedTransportMixtureBuilder builder;

      libMesh::UniquePtr<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::CEACurveFit<libMesh::Real>,
                                                                       Antioch::StatMechThermodynamics<libMesh::Real>,
                                                                       Antioch::SutherlandViscosity<libMesh::Real>,
                                                                       Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                                       Antioch::ConstantLewisDiffusivity<libMesh::Real> > >
        mixture = builder.build_mixture<Antioch::CEACurveFit<libMesh::Real>,
                                        Antioch::StatMechThermodynamics<libMesh::Real>,
                                        Antioch::SutherlandViscosity<libMesh::Real>,
                                        Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                        Antioch::ConstantLewisDiffusivity<libMesh::Real> >( *_input, "TestMaterial" );

      CPPUNIT_ASSERT(mixture);
    }

#ifdef ANTIOCH_HAVE_GSL
    void test_build_cea_kinetic_theory_mix()
    {
      GRINS::AntiochMixtureAveragedTransportMixtureBuilder builder;

      libMesh::UniquePtr<GRINS::AntiochMixtureAveragedTransportMixture<Antioch::CEACurveFit<libMesh::Real>,
                                                                       Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real>,
                                                                       Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                                                       Antioch::KineticsTheoryThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real>,libMesh::Real>,
                                                                       Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> > >
        mixture = builder.build_mixture<Antioch::CEACurveFit<libMesh::Real>,
                                        Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real>,
                                        Antioch::KineticsTheoryViscosity<libMesh::Real,Antioch::GSLSpliner>,
                                        Antioch::KineticsTheoryThermalConductivity<Antioch::IdealGasMicroThermo<Antioch::NASAEvaluator<libMesh::Real,Antioch::CEACurveFit<libMesh::Real> >, libMesh::Real>,libMesh::Real>,
                                        Antioch::MolecularBinaryDiffusion<libMesh::Real,Antioch::GSLSpliner> >( *_input, "TestMaterial" );

      CPPUNIT_ASSERT(mixture);
    }
#endif // ANTIOCH_HAVE_GSL

  private:

    libMesh::UniquePtr<GetPot> _input;
  };


  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochMixtureBuilderBaseTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochConstantTransportMixtureBuilderTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochMixtureAveragedTransportMixtureBuilderTest );

} // end namespace GRINSTesting

#endif // GRINS_HAVE_ANTIOCH
#endif // GRINS_HAVE_CPPUNIT
