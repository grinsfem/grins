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

#include "grins_test_paths.h"
#include "nasa_thermo_test_base.h"
#include "kinetics_test_base.h"
#include "antioch_test_base.h"
#include "cantera_test_base.h"

// GRINS
#include "grins/antioch_evaluator.h"
#include "grins/cantera_evaluator.h"

namespace GRINSTesting
{
  class AirTestBase
  {
  public:

    void init_air_test( unsigned int& N2_idx, unsigned int& O2_idx,
                        unsigned int& N_idx, unsigned int& O_idx,
                        unsigned int& NO_idx,
                        std::vector<unsigned int>& active_species )
    {
      N2_idx = 0;
      O2_idx = 1;
      O_idx  = 4;
      N_idx  = 3;
      NO_idx = 2;

      active_species.resize(5);
      active_species[0] = N2_idx;
      active_species[1] = O2_idx;
      active_species[2] = NO_idx;
      active_species[3] = N_idx;
      active_species[4] = O_idx;
    }

  };

  class AirNASA9Thermo : public AirTestBase,
                         public NASA9ThermoTestBase
  {
  public:

    AirNASA9Thermo()
    {
      this->init_air_test( _N2_idx, _O2_idx, _N_idx, _O_idx, _NO_idx,
                           _active_species );
    }

  };

  class AirKineticsTestBase : public KineticsTestBase,
                              public AirTestBase
  {
  public:

    void init_air_kinetics()
    {
      this->init_air_test( _N2_idx, _O2_idx, _N_idx, _O_idx, _NO_idx,
                           _active_species );

      _n_reactions = 5;

      unsigned int n_species = _active_species.size();

      // All in cal/mol
      _Ea_coeffs.resize(_n_reactions);
      _Ea_coeffs[0] = 224801.3;
      _Ea_coeffs[1] = 117881.7;
      _Ea_coeffs[2] = 149943.0;
      _Ea_coeffs[3] = 85269.6;
      _Ea_coeffs[4] = 38526.0;

      // Convert to J/kmol (since this is what GRINS::Constants::R_universal is in)
      for( unsigned int r = 0; r < _n_reactions; r++ )
        _Ea_coeffs[r] *= 1000.0*4.1868;

      // m^3/kmol
      _preexp_coeffs.resize(_n_reactions);
      _preexp_coeffs[0] = 7.e+18;
      _preexp_coeffs[1] = 2.e+18;
      _preexp_coeffs[2] = 5.e+12;
      _preexp_coeffs[3] = 5.7e+9;
      _preexp_coeffs[4] = 8.4e+09;

      _temp_exp_coeffs.resize(_n_reactions);
      _temp_exp_coeffs[0] = -1.6;
      _temp_exp_coeffs[1] = -1.5;
      _temp_exp_coeffs[2] = 0.0;
      _temp_exp_coeffs[3] = 0.42;
      _temp_exp_coeffs[4] = 0.0;

      _three_body_coeffs.resize(_n_reactions);
      _is_three_body_rxn.resize(_n_reactions,false);
      _is_three_body_rxn[0] = true;
      _three_body_coeffs[0].resize(n_species, 1.0);
      _three_body_coeffs[0][_N_idx] = 4.2857;
      _three_body_coeffs[0][_O_idx] = 4.2857;

      _is_three_body_rxn[1] = true;
      _three_body_coeffs[1].resize(n_species, 1.0);
      _three_body_coeffs[1][_N_idx] = 5.0;
      _three_body_coeffs[1][_O_idx] = 5.0;

      _is_three_body_rxn[2] = true;
      _three_body_coeffs[2].resize(n_species, 1.0);
      _three_body_coeffs[2][_N_idx]  = 22.0;
      _three_body_coeffs[2][_O_idx]  = 22.0;
      _three_body_coeffs[2][_NO_idx] = 22.0;

      _reactant_stoich_coeffs.resize(_n_reactions);
      _product_stoich_coeffs.resize(_n_reactions);

      // N_2 + M <=> 2N + M
      _reactant_stoich_coeffs[0].resize(n_species, 0.0);
      _product_stoich_coeffs[0].resize(n_species, 0.0);
      _reactant_stoich_coeffs[0][_N2_idx] = 1.0;
      _product_stoich_coeffs[0][_N_idx] = 2.0;

      // O_2 + M <=> 2O + M
      _reactant_stoich_coeffs[1].resize(n_species, 0.0);
      _product_stoich_coeffs[1].resize(n_species, 0.0);
      _reactant_stoich_coeffs[1][_O2_idx] = 1.0;
      _product_stoich_coeffs[1][_O_idx] = 2.0;

      // NO + M <=> N + O + M
      _reactant_stoich_coeffs[2].resize(n_species, 0.0);
      _product_stoich_coeffs[2].resize(n_species, 0.0);
      _reactant_stoich_coeffs[2][_NO_idx] = 1.0;
      _product_stoich_coeffs[2][_O_idx] = 1.0;
      _product_stoich_coeffs[2][_N_idx] = 1.0;

      // N2 + 0 <=> NO + N
      _reactant_stoich_coeffs[3].resize(n_species, 0.0);
      _product_stoich_coeffs[3].resize(n_species, 0.0);
      _reactant_stoich_coeffs[3][_N2_idx] = 1.0;
      _reactant_stoich_coeffs[3][_O_idx] = 1.0;
      _product_stoich_coeffs[3][_NO_idx] = 1.0;
      _product_stoich_coeffs[3][_N_idx] = 1.0;

      // NO + 0 <=> O2 + N
      _reactant_stoich_coeffs[4].resize(n_species, 0.0);
      _product_stoich_coeffs[4].resize(n_species, 0.0);
      _reactant_stoich_coeffs[4][_NO_idx] = 1.0;
      _reactant_stoich_coeffs[4][_O_idx] = 1.0;
      _product_stoich_coeffs[4][_O2_idx] = 1.0;
      _product_stoich_coeffs[4][_N_idx] = 1.0;
    }

  };


#ifdef GRINS_HAVE_ANTIOCH

  class AntiochAirNASA9ThermoTest : public CppUnit::TestCase,
                                    public AntiochTestBase,
                                    public AirTestBase,
                                    public NASA9ThermoTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( AntiochAirNASA9ThermoTest );

    CPPUNIT_TEST( test_cp );
    CPPUNIT_TEST( test_hs );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/antioch.in";
      this->init_antioch(input_file, "TestMaterial");

      this->init_air_test( _N2_idx, _O2_idx, _N_idx, _O_idx, _NO_idx,
                           _active_species );

      //this->check_indices(*_antioch_mixture);
    }

    void test_cp()
    {
      std::vector<libMesh::Real> Y(5);
      Y[_N2_idx] = 0.15;
      Y[_O2_idx] = 0.35;
      Y[_NO_idx] = 0.25;
      Y[_O_idx] = 0.2;
      Y[_N_idx] = 0.05;

      this->test_cp_common<GRINS::AntiochMixture,GRINS::AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> > >
        ( *_antioch_mixture, Y, TestingUtils::epsilon()*100 );
    }

    void test_hs()
    {
      this->test_h_common<GRINS::AntiochMixture,GRINS::AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> > >
        ( *_antioch_mixture, TestingUtils::epsilon()*100 );
    }
  };

  class AntiochAirNASA9KineticsTest : public CppUnit::TestCase,
                                      public AntiochTestBase,
                                      public AirKineticsTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( AntiochAirNASA9KineticsTest );

    CPPUNIT_TEST( test_omega_dot );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/antioch.in";
      this->init_antioch(input_file, "TestMaterial");

      this->init_air_kinetics();

      //this->check_indices(*_antioch_mixture);
    }

    void test_omega_dot()
    {
      std::vector<libMesh::Real> Y(5);
      Y[_N2_idx] = 0.15;
      Y[_O2_idx] = 0.35;
      Y[_NO_idx] = 0.25;
      Y[_O_idx] = 0.2;
      Y[_N_idx] = 0.05;

      AirNASA9Thermo thermo;

      this->test_omega_dot_common<GRINS::AntiochMixture,GRINS::AntiochEvaluator<Antioch::CEAEvaluator<libMesh::Real> > >
        ( *_antioch_mixture, thermo, Y, TestingUtils::epsilon()*1e3 );
    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochAirNASA9ThermoTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( AntiochAirNASA9KineticsTest );

#endif // GRINS_HAVE_ANTIOCH


#ifdef GRINS_HAVE_CANTERA

  class CanteraAirNASA9ThermoTest : public CppUnit::TestCase,
                                    public CanteraTestBase,
                                    public AirTestBase,
                                    public NASA9ThermoTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( CanteraAirNASA9ThermoTest );

    CPPUNIT_TEST( test_cp );
    CPPUNIT_TEST( test_hs );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/cantera_chem_thermo.in";
      this->init_cantera(input_file, "TestMaterial");

      this->init_air_test( _N2_idx, _O2_idx, _N_idx, _O_idx, _NO_idx,
                           _active_species );

      //this->check_indices(*_antioch_mixture);
    }

    void test_cp()
    {
      std::vector<libMesh::Real> Y(5);
      Y[_N2_idx] = 0.15;
      Y[_O2_idx] = 0.35;
      Y[_NO_idx] = 0.25;
      Y[_O_idx] = 0.2;
      Y[_N_idx] = 0.05;

      this->test_cp_common<GRINS::CanteraMixture,GRINS::CanteraEvaluator>
        ( *_cantera_mixture, Y, 1.0e-4 );
    }

    void test_hs()
    {
      this->test_h_common<GRINS::CanteraMixture,GRINS::CanteraEvaluator>
        ( *_cantera_mixture, 1.0e-4 );
    }
  };

  class CanteraAirNASA9KineticsTest : public CppUnit::TestCase,
                                      public CanteraTestBase,
                                      public AirKineticsTestBase
  {
  public:
    CPPUNIT_TEST_SUITE( CanteraAirNASA9KineticsTest );

    CPPUNIT_TEST( test_omega_dot );

    CPPUNIT_TEST_SUITE_END();

  public:

    void setUp()
    {
      std::string input_file = std::string(GRINS_TEST_SRCDIR)+"/input_files/cantera_chem_thermo.in";
      this->init_cantera(input_file, "TestMaterial");

      this->init_air_kinetics();

      //this->check_indices(*_antioch_mixture);
    }

    void test_omega_dot()
    {
      std::cout << "Running Cantera test_omega_dot()" << std::endl;

      std::vector<libMesh::Real> Y(5);
      Y[_N2_idx] = 0.15;
      Y[_O2_idx] = 0.35;
      Y[_NO_idx] = 0.25;
      Y[_O_idx] = 0.2;
      Y[_N_idx] = 0.05;

      AirNASA9Thermo thermo;

      this->test_omega_dot_common<GRINS::CanteraMixture,GRINS::CanteraEvaluator>
        ( *_cantera_mixture, thermo, Y, 3.0e-1 );
    }
  };

  CPPUNIT_TEST_SUITE_REGISTRATION( CanteraAirNASA9ThermoTest );
  CPPUNIT_TEST_SUITE_REGISTRATION( CanteraAirNASA9KineticsTest );

#endif // GRINS_HAVE_CANTERA

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT
