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

#ifndef GRINS_NASA_THERMO_TEST_BASE_H
#define GRINS_NASA_THERMO_TEST_BASE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>

#include "species_test_base.h"
#include "thermochem_test_common.h"
#include "testing_utils.h"

namespace GRINSTesting
{
  class NASAThermoTestBase : public SpeciesTestBase
  {
  public:

    virtual libMesh::Real cp_exact( unsigned int species_idx, libMesh::Real T ) = 0;

    virtual libMesh::Real h_exact( unsigned int species_idx, libMesh::Real T ) = 0;

    virtual libMesh::Real s_R_exact( unsigned int species_idx, libMesh::Real T ) = 0;

    virtual libMesh::Real h_RT_exact( unsigned int species_idx, libMesh::Real T ) = 0;

    template<typename ThermoMixture, typename ThermoEvaluator>
    void test_cp_common( ThermoMixture& mixture, const std::vector<libMesh::Real>& Y, libMesh::Real rel_tol )
    {
      CPPUNIT_ASSERT_EQUAL(_active_species.size(), Y.size() );

      ThermoEvaluator evaluator( mixture );

      libMesh::Real R_mix = mixture.R_mix(Y);
      libMesh::Real rho = 1.0e-3;

      libMesh::Real T = 300;
      while( T <= 1000 )
        {
          libMesh::Real P = rho*R_mix*T;

          libMesh::Real cp_mix_computed =  evaluator.cp( T, P, Y );

          std::vector<libMesh::Real> species_cp(_active_species.size(),0.0);
          for( unsigned int s = 0; s < _active_species.size(); s++ )
            species_cp[s] = this->cp_exact(_active_species[s], T);

          libMesh::Real cp_mix_exact =
            ThermochemTestCommon::compute_mass_frac_mixture_prop( species_cp, Y );

          std::stringstream ss;
          ss << T;
          std::string message = "T = "+ss.str();

          libMesh::Real tol = TestingUtils::abs_tol_from_rel_tol( cp_mix_exact, rel_tol );

          CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE( message, cp_mix_exact, cp_mix_computed, tol );

          T += 100.0;
        }
    }

    template<typename ThermoMixture, typename ThermoEvaluator>
    void test_h_common( ThermoMixture& mixture, libMesh::Real rel_tol )
    {
      ThermoEvaluator evaluator( mixture );

      libMesh::Real tol;

      libMesh::Real T = 300;
      while( T <= 1000 )
        {
          std::stringstream ss;
          ss << T;
          std::string message = "T = "+ss.str();

          for( unsigned int s = 0; s < _active_species.size(); s++ )
            {
              libMesh::Real h_exact = this->h_exact(_active_species[s], T);

              libMesh::Real h_computed = evaluator.h_s( T, _active_species[s] );
              tol = TestingUtils::abs_tol_from_rel_tol( h_exact, rel_tol );

              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE( message,
                                                    h_exact,
                                                    h_computed,
                                                    tol );
            }

          T += 100.0;
        }
    }

  protected:

    const std::vector<libMesh::Real>& nasa_coeffs( unsigned int idx )
    {
      if(idx == _N2_idx)
        return _N2_200_1000_coeffs;

      else if( idx == _O2_idx)
        return _O2_200_1000_coeffs;

      else if( idx == _NO_idx)
        return _NO_200_1000_coeffs;

      else if(idx == _O_idx)
        return _O_200_1000_coeffs;

      else if(idx == _N_idx)
        return _N_200_1000_coeffs;

      else
        CPPUNIT_FAIL("Invalid idx for nasa_coeffs");

      // dummy to avoid warning
      return _N2_200_1000_coeffs;
    }

    std::vector<libMesh::Real> _N2_200_1000_coeffs;
    std::vector<libMesh::Real> _O2_200_1000_coeffs;
    std::vector<libMesh::Real> _O_200_1000_coeffs;
    std::vector<libMesh::Real> _NO_200_1000_coeffs;
    std::vector<libMesh::Real> _N_200_1000_coeffs;

  };

  class NASA7ThermoTestBase : public NASAThermoTestBase
  {
  public:

    NASA7ThermoTestBase()
    {
      this->init_N2_coeffs();
      this->init_O2_coeffs();
      this->init_O_coeffs();
      this->init_NO_coeffs();
      this->init_N_coeffs();
    }

    virtual libMesh::Real cp_exact( unsigned int species_idx, libMesh::Real T )
    {
      const std::vector<libMesh::Real> coeffs = this->nasa_coeffs(species_idx);
      return ThermochemTestCommon::nasa7_cp_R_exact(T,
                                                    coeffs[0],
                                                    coeffs[1],
                                                    coeffs[2],
                                                    coeffs[3],
                                                    coeffs[4])*this->R_species(species_idx);
    }

    virtual libMesh::Real h_exact( unsigned int species_idx, libMesh::Real T )
    {
      const std::vector<libMesh::Real> coeffs = this->nasa_coeffs(species_idx);
      return ThermochemTestCommon::nasa7_h_RT_exact(T,
                                                    coeffs[0],
                                                    coeffs[1],
                                                    coeffs[2],
                                                    coeffs[3],
                                                    coeffs[4],
                                                    coeffs[5])*this->R_species(species_idx)*T;
    }

    virtual libMesh::Real s_R_exact( unsigned int species_idx, libMesh::Real T )
    {
      const std::vector<libMesh::Real> coeffs = this->nasa_coeffs(species_idx);
      return ThermochemTestCommon::nasa7_s_R_exact(T,
                                                   coeffs[0],
                                                   coeffs[1],
                                                   coeffs[2],
                                                   coeffs[3],
                                                   coeffs[4],
                                                   coeffs[6]);
    }

    virtual libMesh::Real h_RT_exact( unsigned int species_idx, libMesh::Real T )
    {
      const std::vector<libMesh::Real> coeffs = this->nasa_coeffs(species_idx);
      return ThermochemTestCommon::nasa7_h_RT_exact(T,
                                                    coeffs[0],
                                                    coeffs[1],
                                                    coeffs[2],
                                                    coeffs[3],
                                                    coeffs[4],
                                                    coeffs[5]);
    }

  private:

    void init_N2_coeffs()
    {
      _N2_200_1000_coeffs.resize(7);

      _N2_200_1000_coeffs[0] =  3.53100528E+00;
      _N2_200_1000_coeffs[1] = -1.23660987E-04;
      _N2_200_1000_coeffs[2] = -5.02999437E-07;
      _N2_200_1000_coeffs[3] =  2.43530612E-09;
      _N2_200_1000_coeffs[4] = -1.40881235E-12;
      _N2_200_1000_coeffs[5] = -1.04697628E+03;
      _N2_200_1000_coeffs[6] =  2.96747468E+00;
    }

    void init_O2_coeffs()
    {
      _O2_200_1000_coeffs.resize(7);

      _O2_200_1000_coeffs[0] =  3.78245636E+00;
      _O2_200_1000_coeffs[1] = -2.99673415E-03;
      _O2_200_1000_coeffs[2] =  9.84730200E-06;
      _O2_200_1000_coeffs[3] = -9.68129508E-09;
      _O2_200_1000_coeffs[4] =  3.24372836E-12;
      _O2_200_1000_coeffs[5] = -1.06394356E+03;
      _O2_200_1000_coeffs[6] =  3.65767573E+00;
    }

    void init_O_coeffs()
    {
      _O_200_1000_coeffs.resize(7);

      _O_200_1000_coeffs[0] =  3.16826710E+00;
      _O_200_1000_coeffs[1] = -3.27931884E-03;
      _O_200_1000_coeffs[2] =  6.64306396E-06;
      _O_200_1000_coeffs[3] = -6.12806624E-09;
      _O_200_1000_coeffs[4] =  2.11265971E-12;
      _O_200_1000_coeffs[5] =  2.91222592E+04;
      _O_200_1000_coeffs[6] =  2.05193346E+00;
    }

    void init_NO_coeffs()
    {
      _NO_200_1000_coeffs.resize(7);

      _NO_200_1000_coeffs[0] =  4.21859896E+00;
      _NO_200_1000_coeffs[1] = -4.63988124E-03;
      _NO_200_1000_coeffs[2] =  1.10443049E-05;
      _NO_200_1000_coeffs[3] = -9.34055507E-09;
      _NO_200_1000_coeffs[4] =  2.80554874E-12;
      _NO_200_1000_coeffs[5] =  9.84509964E+03;
      _NO_200_1000_coeffs[6] =  2.28061001E+00;
    }

    void init_N_coeffs()
    {
      _N_200_1000_coeffs.resize(7);

      _N_200_1000_coeffs[0] =  2.50000000E+00;
      _N_200_1000_coeffs[1] =  0.0;
      _N_200_1000_coeffs[2] =  0.0;
      _N_200_1000_coeffs[3] =  0.0;
      _N_200_1000_coeffs[4] =  0.0;
      _N_200_1000_coeffs[5] =  5.61046378E+04;
      _N_200_1000_coeffs[6] =  4.19390932E+00;
    }

  };

  class NASA9ThermoTestBase : public NASAThermoTestBase
  {
  public:

    NASA9ThermoTestBase()
    {
      this->init_N2_coeffs();
      this->init_O2_coeffs();
      this->init_O_coeffs();
      this->init_NO_coeffs();
      this->init_N_coeffs();
    }

    virtual libMesh::Real cp_exact( unsigned int species_idx, libMesh::Real T )
    {
      const std::vector<libMesh::Real> coeffs = this->nasa_coeffs(species_idx);
      return ThermochemTestCommon::nasa9_cp_R_exact(T,
                                                    coeffs[0],
                                                    coeffs[1],
                                                    coeffs[2],
                                                    coeffs[3],
                                                    coeffs[4],
                                                    coeffs[5],
                                                    coeffs[6])*this->R_species(species_idx);
    }

    virtual libMesh::Real h_exact( unsigned int species_idx, libMesh::Real T )
    {
      const std::vector<libMesh::Real> coeffs = this->nasa_coeffs(species_idx);
      return ThermochemTestCommon::nasa9_h_RT_exact(T,
                                                    coeffs[0],
                                                    coeffs[1],
                                                    coeffs[2],
                                                    coeffs[3],
                                                    coeffs[4],
                                                    coeffs[5],
                                                    coeffs[6],
                                                    coeffs[7])*this->R_species(species_idx)*T;
    }

    virtual libMesh::Real s_R_exact( unsigned int species_idx, libMesh::Real T )
    {
      const std::vector<libMesh::Real> coeffs = this->nasa_coeffs(species_idx);
      return ThermochemTestCommon::nasa9_s_R_exact(T,
                                                   coeffs[0],
                                                   coeffs[1],
                                                   coeffs[2],
                                                   coeffs[3],
                                                   coeffs[4],
                                                   coeffs[5],
                                                   coeffs[6],
                                                   coeffs[8]);
    }

    virtual libMesh::Real h_RT_exact( unsigned int species_idx, libMesh::Real T )
    {
      const std::vector<libMesh::Real> coeffs = this->nasa_coeffs(species_idx);
      return ThermochemTestCommon::nasa9_h_RT_exact(T,
                                                    coeffs[0],
                                                    coeffs[1],
                                                    coeffs[2],
                                                    coeffs[3],
                                                    coeffs[4],
                                                    coeffs[5],
                                                    coeffs[6],
                                                    coeffs[7]);
    }

  private:

    void init_N2_coeffs()
    {
      _N2_200_1000_coeffs.resize(9);

      _N2_200_1000_coeffs[0] =  2.21037122e+04;
      _N2_200_1000_coeffs[1] = -3.81846145e+02;
      _N2_200_1000_coeffs[2] =  6.08273815e+00;
      _N2_200_1000_coeffs[3] = -8.53091381e-03;
      _N2_200_1000_coeffs[4] =  1.38464610e-05;
      _N2_200_1000_coeffs[5] = -9.62579293e-09;
      _N2_200_1000_coeffs[6] =  2.51970560e-12;
      _N2_200_1000_coeffs[7] =  7.10845911e+02;
      _N2_200_1000_coeffs[8] = -1.07600320e+01;
    }

    void init_O2_coeffs()
    {
      _O2_200_1000_coeffs.resize(9);

      _O2_200_1000_coeffs[0] = -3.42556269e+04;
      _O2_200_1000_coeffs[1] =  4.84699986e+02;
      _O2_200_1000_coeffs[2] =  1.11901159e+00;
      _O2_200_1000_coeffs[3] =  4.29388743e-03;
      _O2_200_1000_coeffs[4] = -6.83627313e-07;
      _O2_200_1000_coeffs[5] = -2.02337478e-09;
      _O2_200_1000_coeffs[6] =  1.03904064e-12;
      _O2_200_1000_coeffs[7] = -3.39145434e+03;
      _O2_200_1000_coeffs[8] =  1.84969912e+01;
    }

    void init_O_coeffs()
    {
      _O_200_1000_coeffs.resize(9);

      _O_200_1000_coeffs[0] = -7.95361130e+03;
      _O_200_1000_coeffs[1] =  1.60717779e+02;
      _O_200_1000_coeffs[2] =  1.96622644e+00;
      _O_200_1000_coeffs[3] =  1.01367031e-03;
      _O_200_1000_coeffs[4] = -1.11041542e-06;
      _O_200_1000_coeffs[5] =  6.51750750e-10;
      _O_200_1000_coeffs[6] = -1.58477925e-13;
      _O_200_1000_coeffs[7] =  2.84036244e+04;
      _O_200_1000_coeffs[8] =  8.40424182e+00;
    }

    void init_NO_coeffs()
    {
      _NO_200_1000_coeffs.resize(9);

      _NO_200_1000_coeffs[0] = -1.14391658e+04;
      _NO_200_1000_coeffs[1] =  1.53646774e+02;
      _NO_200_1000_coeffs[2] =  3.43146865e+00;
      _NO_200_1000_coeffs[3] = -2.66859213e-03;
      _NO_200_1000_coeffs[4] =  8.48139877e-06;
      _NO_200_1000_coeffs[5] = -7.68511079e-09;
      _NO_200_1000_coeffs[6] =  2.38679758e-12;
      _NO_200_1000_coeffs[7] =  9.09794974e+03;
      _NO_200_1000_coeffs[8] =  6.72872795e+00;
    }

    void init_N_coeffs()
    {
      _N_200_1000_coeffs.resize(9);

      _N_200_1000_coeffs[0] =  0.0;
      _N_200_1000_coeffs[1] =  0.0;
      _N_200_1000_coeffs[2] =  2.5;
      _N_200_1000_coeffs[3] =  0.0;
      _N_200_1000_coeffs[4] =  0.0;
      _N_200_1000_coeffs[5] =  0.0;
      _N_200_1000_coeffs[6] =  0.0;
      _N_200_1000_coeffs[7] =  5.61046378e+04;
      _N_200_1000_coeffs[8] =  4.19390932e+00;
    }

  };

} // end namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT

#endif // GRINS_NASA_THERMO_TEST_BASE_H
