//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2019 Paul T. Bauman, Roy H. Stogner
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


// GRINS
#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// C++
#include <iomanip>
#include <vector>
#include <fstream>

// libMesh
#include "libmesh/getpot.h"

// GRINS
#include "grins/antioch_mixture_averaged_transport_evaluator.h"
#include "grins/materials_parsing.h"
#include "grins/physics_naming.h"

template<typename KineticsThermo, typename Thermo, typename Viscosity, typename Conductivity, typename Diffusivity>
int do_transport_eval( const GetPot& input )
{
  GRINS::AntiochMixtureAveragedTransportMixture<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity> mixture(input,GRINS::MaterialsParsing::material_name(input,GRINS::PhysicsNaming::reacting_low_mach_navier_stokes()));

  GRINS::AntiochMixtureAveragedTransportEvaluator<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity> evaluator(mixture);

  libMesh::Real T0 = input( "Conditions/T0", 300.0 );
  libMesh::Real T1 = input( "Conditions/T1", 300.0 );
  libMesh::Real T_inc = input( "Conditions/T_increment", 100.0 );

  libMesh::Real rho = input( "Conditions/density", 1.0e-3 );

  const unsigned int n_species = mixture.n_species();

  std::vector<libMesh::Real> Y(n_species);
  if( input.vector_variable_size( "Conditions/mass_fractions" ) != n_species )
    {
      std::cerr << "Error: mass fractions size not consistent with n_species"
                << std::endl;
      libmesh_error();
    }

  for( unsigned int s = 0; s < n_species; s++ )
    {
      Y[s] = input( "Conditions/mass_fractions", 0.0, s );
    }

  libMesh::Real T = T0;

  libMesh::Real p0 = rho*T*evaluator.R_mix(Y);

  std::ofstream output;
  output.open( "transport.dat", std::ios::trunc );

  output << "# Species names" << std::endl;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      output << mixture.species_name( s ) << " ";
    }
  output << std::endl;
  output << "# T [K]    mu         k           D[s]" << std::endl;

  output.close();

  while( T < T1 )
    {
      output.open( "transport.dat", std::ios::app );
      output << std::scientific << std::setprecision(16);
      output << T << " ";

      libMesh::Real mu;
      libMesh::Real k;
      std::vector<libMesh::Real> D(n_species);
      evaluator.mu_and_k_and_D( T, rho, evaluator.cp(T,p0,Y), Y, mu, k, D );

      output <<  mu << " ";
      output << k << " ";

      for( unsigned int s = 0; s< n_species; s++ )
        {
          output << D[s] << " ";
        }
      output << std::endl;
      output.close();

      T += T_inc;
    }

  return 0;
}

int main(int argc, char* argv[])
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input file." << std::endl;
      exit(1);
    }

  GetPot input( argv[1] );

  std::string mixing_model = input( "Physics/Antioch/mixing_model", "wilke");
  std::string thermo_model = input( "Physics/Antioch/thermo_model", "stat_mech");
  std::string viscosity_model = input( "Physics/Antioch/viscosity_model", "sutherland");
  std::string conductivity_model = input( "Physics/Antioch/conductivity_model", "eucken");
  std::string diffusivity_model = input( "Physics/Antioch/diffusivity_model", "constant_lewis");
  std::string kinetics_thermo_model = input( "Physics/Antioch/kinetics_thermo_model", "cea");

  int return_flag = 0;

  if( mixing_model == std::string("wilke") )
    {
      if( kinetics_thermo_model == std::string("cea") )
        {
          if( thermo_model == std::string("stat_mech") )
            {
              if( diffusivity_model == std::string("constant_lewis") )
                {
                  if( conductivity_model == std::string("eucken") )
                    {
                      if( viscosity_model == std::string("sutherland") )
                        {
                          return_flag = do_transport_eval<Antioch::CEACurveFit<libMesh::Real>,
                                                          Antioch::StatMechThermodynamics<libMesh::Real>,
                                                          Antioch::SutherlandViscosity<libMesh::Real>,
                                                          Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                          Antioch::ConstantLewisDiffusivity<libMesh::Real> >(input);
                        }
                      else if( viscosity_model == std::string("blottner") )
                        {
                          return_flag = do_transport_eval<Antioch::CEACurveFit<libMesh::Real>,
                                                          Antioch::StatMechThermodynamics<libMesh::Real>,
                                                          Antioch::BlottnerViscosity<libMesh::Real>,
                                                          Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                                                          Antioch::ConstantLewisDiffusivity<libMesh::Real> >(input);
                        }
                      else
                        {
                          std::cerr << "Error: Unknown viscosity_model "
                                    << viscosity_model << "!"  << std::endl;
                          return_flag = 1;
                        }
                    }
                  else
                    {
                      std::cerr << "Error: Unknown conductivity_model "
                                << conductivity_model << "!"  << std::endl;
                      return_flag = 1;
                    }
                }
              else
                {
                  std::cerr << "Error: Unknown diffusivity_model "
                            << diffusivity_model << "!"  << std::endl;
                  return_flag = 1;
                }
            }
          else
            {
              std::cerr << "Error: Unknown thermo_model "
                        << thermo_model << "!"  << std::endl;
              return_flag = 1;
            }
        }
      else
        {
          std::cerr << "Error: Unknown kinetics_thermo_model "
                    << thermo_model << "!"  << std::endl;
          return_flag = 1;
        }
    }
  else
    {
      std::cerr << "Error: Unknown mixing_model "
                << mixing_model << "!" << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

#endif //GRINS_HAVE_ANTIOCH
