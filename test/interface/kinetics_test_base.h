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

#ifndef GRINS_KINETICS_TEST_BASE_H
#define GRINS_KINETICS_TEST_BASE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_CPPUNIT

#include <cppunit/extensions/HelperMacros.h>

#include "species_test_base.h"
#include "thermochem_test_common.h"
#include "testing_utils.h"

#include "libmesh/libmesh_common.h"

namespace GRINSTesting
{
  class KineticsTestBase : public SpeciesTestBase
  {
  public:

    template<typename ThermoMixture, typename ThermoEvaluator>
    void test_omega_dot_common( ThermoMixture& mixture, NASAThermoTestBase& thermo_funcs,
                                const std::vector<libMesh::Real>& Y, libMesh::Real rel_tol )
    {
      CPPUNIT_ASSERT_EQUAL(_active_species.size(), Y.size() );

      ThermoEvaluator evaluator( mixture );

      libMesh::Real rho = 1.0e-3;

      std::vector<libMesh::Real> molar_densities(_active_species.size(),0.0);
      for( unsigned int s = 0; s < _active_species.size(); s++ )
        molar_densities[s] = rho*Y[s]/this->molar_mass(s);

      libMesh::Real T = 300;
      while( T <= 1000 )
        {
          std::vector<libMesh::Real> forward_rates, backward_rates;
          this->compute_reaction_rates( T, molar_densities, thermo_funcs, forward_rates, backward_rates );

          std::vector<libMesh::Real> omega_dot_exact(_active_species.size(), 0.0);
          std::vector<libMesh::Real> omega_dot_computed(_active_species.size());

          evaluator.omega_dot( T, rho, Y, omega_dot_computed );

          for( unsigned int s = 0; s < _active_species.size(); s++ )
            for( unsigned int r = 0; r < _n_reactions; r++ )
              {
                unsigned int species = _active_species[s];
                libMesh::Real stoich = _product_stoich_coeffs[r][species] - _reactant_stoich_coeffs[r][species];

                omega_dot_exact[s] += this->molar_mass(species)*stoich*( forward_rates[r] - backward_rates[r] );
              }

          for( unsigned int s = 0; s < _active_species.size(); s++ )
            {
              unsigned int species = _active_species[s];

              std::stringstream ss;
              ss << T;
              std::string message = "T = "+ss.str();

              ss.str(std::string());
              ss << s;
              message += ", species = "+mixture.species_name(species);

              /*
              std::cout << message
                        << ", omega_dot_exact = " << omega_dot_exact[s]
                        << ", omega_dot_computed = " << omega_dot_computed[s]
                        << std::endl;
              */

              libMesh::Real tol = TestingUtils::abs_tol_from_rel_tol( omega_dot_exact[s], rel_tol );
              CPPUNIT_ASSERT_DOUBLES_EQUAL_MESSAGE( message,
                                                    omega_dot_exact[s],
                                                    omega_dot_computed[s],
                                                    tol );
            }

          T += 100.0;
        }
    }

  protected:

    unsigned int _n_reactions;

    std::vector<libMesh::Real> _Ea_coeffs;
    std::vector<libMesh::Real> _preexp_coeffs;
    std::vector<libMesh::Real> _temp_exp_coeffs;

    std::vector<std::vector<libMesh::Real> > _three_body_coeffs;
    std::vector<bool> _is_three_body_rxn;

    std::vector<std::vector<libMesh::Real> > _reactant_stoich_coeffs;
    std::vector<std::vector<libMesh::Real> > _product_stoich_coeffs;

  private:

    void compute_reaction_rates( libMesh::Real T,
                                 const std::vector<libMesh::Real>& molar_densities,
                                 NASAThermoTestBase& thermo_funcs,
                                 std::vector<libMesh::Real>& forward_rates,
                                 std::vector<libMesh::Real>& backward_rates )
    {
      forward_rates.resize(_n_reactions,1.0);
      backward_rates.resize(_n_reactions,1.0);

      for( unsigned int r = 0; r < _n_reactions; r++ )
        {
          libMesh::Real fwd_rate_coeff = ThermochemTestCommon::arrhenius_rate( _preexp_coeffs[r],
                                                                               _temp_exp_coeffs[r],
                                                                               _Ea_coeffs[r],
                                                                               T );

          libMesh::Real Keq = this->eq_constant(T,
                                                _reactant_stoich_coeffs[r],
                                                _product_stoich_coeffs[r],
                                                thermo_funcs);

          libMesh::Real bkwd_rate_coeff = fwd_rate_coeff/Keq;

          forward_rates[r] *= fwd_rate_coeff;
          backward_rates[r] *= bkwd_rate_coeff;

          for( unsigned int s = 0; s < _active_species.size(); s++ )
            {
              forward_rates[r] *= std::pow(molar_densities[s],_reactant_stoich_coeffs[r][s]);

              backward_rates[r] *= std::pow(molar_densities[s],_product_stoich_coeffs[r][s]);
            }

          if( _is_three_body_rxn[r] )
            {
              libMesh::Real M = ThermochemTestCommon::compute_third_body_molar_density( molar_densities,
                                                                                        _three_body_coeffs[r] );
              forward_rates[r] *= M;
              backward_rates[r] *= M;
            }
        }
    }

    libMesh::Real eq_constant( libMesh::Real T,
                               std::vector<libMesh::Real>& reactant_stoich_coeffs,
                               std::vector<libMesh::Real>& product_stoich_coeffs,
                               NASAThermoTestBase& thermo_funcs )
    {
      CPPUNIT_ASSERT_EQUAL( reactant_stoich_coeffs.size(), product_stoich_coeffs.size() );

      libMesh::Real coeff_sum = 0.0;
      libMesh::Real exp_sum = 0.0;

      for( unsigned int s = 0; s < reactant_stoich_coeffs.size(); s++ )
        {
          libMesh::Real stoich = (product_stoich_coeffs[s] - reactant_stoich_coeffs[s]);

          coeff_sum += stoich;
          exp_sum += stoich*( thermo_funcs.h_RT_exact(s,T) - thermo_funcs.s_R_exact(s,T) );
        }

      libMesh::Real P_RT = 1e5/(GRINS::Constants::R_universal*T);

      return std::pow(P_RT,coeff_sum)*std::exp(-exp_sum);
    }

  };

} // namespace GRINSTesting

#endif // GRINS_HAVE_CPPUNIT

#endif // GRINS_KINETICS_TEST_BASE_H
