//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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

#ifdef GRINS_HAVE_CANTERA

// This class
#include "grins/cantera_kinetics.h"

// GRINS
#include "grins/cached_values.h"
#include "grins/cantera_mixture.h"

// libMesh
#include "libmesh/getpot.h"

// Cantera (with compiler warnings disabled)
#include "libmesh/ignore_warnings.h"
#include "cantera/IdealGasMix.h"
#include "libmesh/restore_warnings.h"

namespace GRINS
{

  CanteraKinetics::CanteraKinetics( CanteraMixture& mixture )
    :  _cantera_gas( mixture.get_chemistry() )
  {
    return;
  }

  CanteraKinetics::~CanteraKinetics()
  {
    return;
  }

  void CanteraKinetics::omega_dot( const CachedValues& cache,
				   unsigned int qp,
				   std::vector<libMesh::Real>& omega_dot ) const
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const libMesh::Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    this->omega_dot_TPY( T, P, Y, omega_dot );

    return;
  }

  void CanteraKinetics::omega_dot_TPY( const libMesh::Real T, const libMesh::Real P,
                                       const std::vector<libMesh::Real>& mass_fractions,
                                       std::vector<libMesh::Real>& omega_dot ) const
  {
    libmesh_assert_equal_to( mass_fractions.size(), omega_dot.size() );
    libmesh_assert_equal_to( mass_fractions.size(), _cantera_gas.nSpecies() );
    libmesh_assert_greater(T,0.0);
    libmesh_assert_greater(P,0.0);
    
    {
      /*! \todo Need to make sure this will work in a threaded environment.
	Not sure if we will get thread lock here or not. */
      libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);
      
      try
	{
	  _cantera_gas.setState_TPY(T, P, &mass_fractions[0]);
	  _cantera_gas.getNetProductionRates(&omega_dot[0]);
	}
      catch(Cantera::CanteraError)
	{
	  Cantera::showErrors(std::cerr);
	  libmesh_error();
	}
    }

#ifdef DEBUG
    for( unsigned int s = 0; s < omega_dot.size(); s++ )
      {
        if( libMesh::libmesh_isnan(omega_dot[s]) )
          {
            std::cout << "T = " << T << std::endl
                      << "P = " << P << std::endl;
            for( unsigned int s = 0; s < omega_dot.size(); s++ )
              {
                std::cout << "Y[" << s << "] = " << mass_fractions[s] << std::endl;
              }
            for( unsigned int s = 0; s < omega_dot.size(); s++ )
              {
                std::cout << "omega_dot[" << s << "] = " << omega_dot[s] << std::endl;
              }

            libmesh_error();
          }
      }
#endif

    for( unsigned int s = 0; s < omega_dot.size(); s++ )
      {
        // convert [kmol/m^3-s] to [kg/m^3-s]
        omega_dot[s] *= this->_cantera_gas.molecularWeight(s);
      }

    return;
  }

  void CanteraKinetics::omega_dot_TRY( const libMesh::Real& T, const libMesh::Real rho,
                                       const std::vector<libMesh::Real>& mass_fractions,
                                       std::vector<libMesh::Real>& omega_dot ) const
  {
    libmesh_assert_equal_to( mass_fractions.size(), omega_dot.size() );
    libmesh_assert_equal_to( mass_fractions.size(), _cantera_gas.nSpecies() );
    libmesh_assert_greater(T,0.0);
    libmesh_assert_greater(rho,0.0);

    {
      /*! \todo Need to make sure this will work in a threaded environment.
	Not sure if we will get thread lock here or not. */
      libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

      try
	{
	  _cantera_gas.setState_TRY(T, rho, &mass_fractions[0]);
	  _cantera_gas.getNetProductionRates(&omega_dot[0]);
	}
      catch(Cantera::CanteraError)
	{
	  Cantera::showErrors(std::cerr);
	  libmesh_error();
	}
    }

    for( unsigned int s = 0; s < omega_dot.size(); s++ )
      {
        // convert [kmol/m^3-s] to [kg/m^3-s]
        omega_dot[s] *= this->_cantera_gas.molecularWeight(s);
      }

    return;
  }

} // namespace GRINS

#endif //GRINS_HAVE_CANTERA
