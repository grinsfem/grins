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


#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// C++
#include <iomanip>
#include <limits>
#include <vector>

// GRINS
#include "grins/antioch_mixture_averaged_transport_mixture.h"
#include "grins/antioch_mixture_averaged_transport_evaluator.h"
#include "grins/cached_values.h"
#include "grins/materials_parsing.h"
#include "grins/physics_naming.h"
#include "grins/antioch_mixture_averaged_transport_mixture_builder.h"

// libMesh
#include "libmesh/getpot.h"

int test_generic( const libMesh::Real value, const libMesh::Real value_reg, const std::string& name )
{
  int return_flag = 0;

  const double tol = std::numeric_limits<double>::epsilon()*10;

  const double rel_error = std::fabs( (value - value_reg)/value_reg );

  if( rel_error > tol )
    {
      return_flag = 1;
      std::cout << "Mismatch in "+name  << std::endl
                << name+" = " << value << std::endl
                << name+"_reg = " << value_reg << std::endl
                << "rel_error = " << rel_error << std::endl;
    }

  return return_flag;
}


template<typename Viscosity>
int test_mu( const libMesh::Real mu );

template<typename Thermo, typename Conductivity>
int test_k( const libMesh::Real k );

template<typename Thermo, typename Viscosity, typename Conductivity, typename Diffusivity>
int test_D( const std::vector<libMesh::Real>& D );

template<>
int test_mu<Antioch::BlottnerViscosity<libMesh::Real> >( const libMesh::Real mu )
{
  double mu_reg = 4.5123309407810213e-05;

  return test_generic(mu,mu_reg,"mu");
}

template<>
int test_k<Antioch::StatMechThermodynamics<libMesh::Real>,
           Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> > >( const libMesh::Real k )
{
  double k_reg = 8.0102737519532258e-02;

  return test_generic(k,k_reg,"k");
}

template<>
int test_D<Antioch::StatMechThermodynamics<libMesh::Real>,
           Antioch::BlottnerViscosity<libMesh::Real>,
           Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
           Antioch::ConstantLewisDiffusivity<libMesh::Real> >( const std::vector<libMesh::Real>& D )
{
  std::vector<libMesh::Real> D_reg(5);
  D_reg[0] = 4.6482311273552429e-02;
  D_reg[1] = D_reg[0];
  D_reg[2] = D_reg[0];
  D_reg[3] = D_reg[0];
  D_reg[4] = D_reg[0];

  int return_flag_tmp = 0;
  int return_flag = 0;
  for( unsigned int s = 0; s < 5; s++ )
    {
      return_flag_tmp = test_generic(D[s],D_reg[s],"D_s");
      if( return_flag_tmp != 0 ) return_flag = 1;
    }

  return return_flag;
}

template<typename KineticsThermo, typename Thermo, typename Viscosity, typename Conductivity, typename Diffusivity>
int test_evaluator( const GetPot& input )
{
  GRINS::AntiochMixtureAveragedTransportMixtureBuilder builder;

  std::unique_ptr<GRINS::AntiochMixtureAveragedTransportMixture<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity> >
    mixture_ptr = builder.build_mixture<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity>
    (input,"TestMaterial");


  const GRINS::AntiochMixtureAveragedTransportMixture<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity> &
    mixture = *mixture_ptr;

  GRINS::AntiochMixtureAveragedTransportEvaluator<KineticsThermo,Thermo,Viscosity,Conductivity,Diffusivity> evaluator(mixture);

  const libMesh::Real T = 1000;

  const libMesh::Real rho = 1.0e-3;

  const unsigned int n_species = 5;

  std::vector<libMesh::Real> Y(n_species,0.2);

  libMesh::Real p0 = rho*T*evaluator.R_mix(Y);
  libMesh::Real mu = 0.0;
  libMesh::Real k = 0.0;
  std::vector<libMesh::Real> D(n_species,0.0);

  evaluator.mu_and_k_and_D( T, rho, evaluator.cp(T,p0,Y), Y, mu, k, D );

  std::cout << std::scientific << std::setprecision(16)
            << "mu = " << mu << std::endl;

  std::cout << std::scientific << std::setprecision(16)
            << "k = " << k << std::endl;

  for( unsigned int i = 0; i < n_species; i++ )
    {
      std::cout << std::scientific << std::setprecision(16)
                << "D(" << mixture.species_name(i) << ") = " << D [i] << std::endl;
    }

  int return_flag = 0;

  int return_flag_temp = 0;

  return_flag_temp = test_mu<Viscosity>( mu );
  if( return_flag_temp != 0 ) return_flag = 1;

  return_flag_temp = test_k<Thermo,Conductivity>( k );
  if( return_flag_temp != 0 ) return_flag = 1;

  return_flag_temp = test_D<Thermo,Viscosity,Conductivity,Diffusivity>( D );
  if( return_flag_temp != 0 ) return_flag = 1;

  return return_flag;
}


int main( int argc, char* argv[] )
{
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify input file." << std::endl;
      exit(1);
    }

  GetPot input( argv[1] );

  int return_flag = 0;

  std::cout << std::endl <<  "Running StatMesh, Blottner, Eucken, Constant Lewis regression test." << std::endl;
  return_flag = test_evaluator<Antioch::CEACurveFit<libMesh::Real>,
                               Antioch::StatMechThermodynamics<libMesh::Real>,
                               Antioch::BlottnerViscosity<libMesh::Real>,
                               Antioch::EuckenThermalConductivity<Antioch::StatMechThermodynamics<libMesh::Real> >,
                               Antioch::ConstantLewisDiffusivity<libMesh::Real> >(input);

  return return_flag;
}

#else //GRINS_HAVE_ANTIOCH
int main()
{
  // automake expects 77 for a skipped test
  return 77;
}
#endif
