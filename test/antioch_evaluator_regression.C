//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
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
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// C++
#include <iomanip>
#include <limits>
#include <vector>

// GRINS
#include "grins/antioch_mixture.h"
#include "grins/antioch_evaluator.h"
#include "grins/cached_values.h"

// libMesh
#include "libmesh/getpot.h"

// Antioch
#include "antioch/cea_evaluator.h"
#include "antioch/stat_mech_thermo.h"

// Boost
#include "boost/math/special_functions/fpclassify.hpp" //isnan

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

template<typename Thermo>
int test_cp( const libMesh::Real cp );

template<typename Thermo>
int test_cv( const libMesh::Real cv );

template<typename Thermo>
int test_h_s( const libMesh::Real h_s );

template<>
int test_cp<Antioch::CEAEvaluator<libMesh::Real> >( const libMesh::Real cp )
{
  double cp_reg = 1.2361869971209990e+03;

  return test_generic(cp,cp_reg,"cp");
}

template<>
int test_cv<Antioch::CEAEvaluator<libMesh::Real> >( const libMesh::Real cv )
{
  double cv_reg = 8.4681056933423179e+02;

  return test_generic(cv,cv_reg,"cv");
}

template<>
int test_h_s<Antioch::CEAEvaluator<libMesh::Real> >( const libMesh::Real h_s )
{
  double h_s_reg = 7.6606764036494098e+05;

  return test_generic(h_s,h_s_reg,"h_s");
}

template<>
int test_cp<Antioch::StatMechThermodynamics<libMesh::Real> >( const libMesh::Real cp )
{
  double cp_reg = 1.2309250447693457e+03;

  return test_generic(cp,cp_reg,"cp");
}

template<>
int test_cv<Antioch::StatMechThermodynamics<libMesh::Real> >( const libMesh::Real cv )
{
  double cv_reg = 8.4154861698257866e+02;

  return test_generic(cv,cv_reg,"cv");
}

template<>
int test_h_s<Antioch::StatMechThermodynamics<libMesh::Real> >( const libMesh::Real h_s )
{
  double h_s_reg = 7.6397682090011111e+05;

  return test_generic(h_s,h_s_reg,"h_s");
}


template <typename Thermo>
int test_evaluator( const GRINS::AntiochMixture& antioch_mixture )
{
  GRINS::AntiochEvaluator<Thermo> antioch_evaluator( antioch_mixture );

  libMesh::Real T = 1000;

  libMesh::Real rho = 1.0e-3;

  const unsigned int n_species = 5;

  std::vector<libMesh::Real> Y(n_species,0.2);

  libMesh::Real R_mix = antioch_mixture.R_mix(Y);

  GRINS::CachedValues cache;

  cache.add_quantity(GRINS::Cache::TEMPERATURE);
  std::vector<double> Tqp(1,T);
  cache.set_values(GRINS::Cache::TEMPERATURE, Tqp);

  cache.add_quantity(GRINS::Cache::MIXTURE_DENSITY);
  std::vector<double> rhoqp(1,rho);
  cache.set_values(GRINS::Cache::MIXTURE_DENSITY, rhoqp);

  cache.add_quantity(GRINS::Cache::MIXTURE_GAS_CONSTANT);
  std::vector<double> Rqp(1,R_mix);
  cache.set_values(GRINS::Cache::MIXTURE_GAS_CONSTANT, Rqp);

  cache.add_quantity(GRINS::Cache::MASS_FRACTIONS);
  std::vector<std::vector<double> > Yqp(1,Y);
  cache.set_vector_values(GRINS::Cache::MASS_FRACTIONS, Yqp);

  std::vector<libMesh::Real> omega_dot(n_species,0.0);

  libMesh::Real cp =  antioch_evaluator.cp( cache, 0 );

  std::cout << std::scientific << std::setprecision(16)
            << "cp = " << cp << std::endl;

  libMesh::Real cv =  antioch_evaluator.cv( cache, 0 );

  std::cout << std::scientific << std::setprecision(16)
            << "cv = " << cv << std::endl;

  libMesh::Real h_s =  antioch_evaluator.h_s( cache, 0, 0 );

  std::cout << std::scientific << std::setprecision(16)
            << "h_s = " << h_s << std::endl;

  int return_flag = 0;
  int return_flag_temp = 0;

  return_flag_temp = test_cp<Thermo>( cp );
  if( return_flag_temp != 0 ) return_flag = 1;

  return_flag_temp = test_cv<Thermo>( cv );
  if( return_flag_temp != 0 ) return_flag = 1;
  
  return_flag_temp = test_h_s<Thermo>( h_s );
  if( return_flag_temp != 0 ) return_flag = 1;

  antioch_evaluator.omega_dot( cache, 0, omega_dot );

  for( unsigned int i = 0; i < n_species; i++ )
    {
      std::cout << std::scientific << std::setprecision(16) 
                << "omega_dot(" << i << ") = " << omega_dot[i] << std::endl;
    }

  
  double tol = std::numeric_limits<double>::epsilon()*10;

  // Check that omega_dot sums to 1
  libMesh::Real sum = 0.0;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      sum += omega_dot[s];
    }

  if( std::fabs(sum) > tol )
    {
      std::cout << "Error: Sum of mass sources not equal to zero!" << std::endl
                << "sum = " << sum << std::endl;
      return_flag = 1;
    }

  bool omega_is_nan = false;
  for( unsigned int s = 0; s < n_species; s++ )
    {
      if( boost::math::isnan(omega_dot[s]) )
        {
          omega_is_nan = true;
          return_flag = 1;
        }
    }

  if( omega_is_nan )
    {
      std::cerr << "Error: Detected NAN's for omega_dot!" << std::endl;
    }

  // Now check against regression values
  std::vector<libMesh::Real> omega_dot_reg(n_species,0.0);
  omega_dot_reg[0] =  4.3563269373325170e-01;
  omega_dot_reg[1] = -3.6798048754802180e+00;
  omega_dot_reg[2] =  2.9971144322942522e+00;
  omega_dot_reg[3] = -1.8347122381073480e+00;
  omega_dot_reg[4] =  2.0817699875600622e+00;

  for( unsigned int s = 0; s < n_species; s++ )
    {
      if( std::fabs( (omega_dot[s] - omega_dot_reg[s])/omega_dot_reg[s] ) > tol )
	{
	  std::cerr << "Error: Mismatch in omega_dot." << std::endl
		    << std::setprecision(16) << std::scientific
		    << "s = " << s << std::endl
		    << "omega_dot = " << omega_dot[s] << std::endl
		    << "omega_dot_reg = " << omega_dot_reg[s] << std::endl;
	  return_flag = 1;
	}
    }

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
  
  GRINS::AntiochMixture antioch_mixture(input);

  int return_flag = 0;

  std::cout << "Running AntiochCEAThermo regression test." << std::endl;
  return_flag = test_evaluator<Antioch::CEAEvaluator<libMesh::Real> >( antioch_mixture );

  std::cout << std::endl <<  "Running AntiochStatMechThermo regression test." << std::endl;
  return_flag = test_evaluator<Antioch::StatMechThermodynamics<libMesh::Real> >( antioch_mixture );

  return return_flag;
}

#else //GRINS_HAVE_ANTIOCH
int main()
{
  // automake expects 77 for a skipped test
  return 77;
}
#endif
