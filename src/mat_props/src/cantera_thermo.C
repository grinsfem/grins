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

#include "grins/cantera_thermo.h"

#ifdef GRINS_HAVE_CANTERA
namespace
{
  libMesh::Threads::spin_mutex thermo_mutex;
}

namespace GRINS
{

  CanteraThermodynamics::CanteraThermodynamics( const GetPot& input, 
						const ChemicalMixture& chem_mixture )
    : _chem_mixture(chem_mixture),
      _cantera_gas( CanteraSingleton::cantera_instance(input) )
  {
    return;
  }

  CanteraThermodynamics::~CanteraThermodynamics()
  {
    return;
  }

  libMesh::Real CanteraThermodynamics::cp( const CachedValues& cache, unsigned int qp ) const
  {
    libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const libMesh::Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];
    
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    libMesh::Real cp = 0.0;

    {
      //libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

      /*! \todo Need to make sure this will work in a threaded environment.
	Not sure if we will get thread lock here or not. */
      try
	{
	  _cantera_gas.setState_TPY( T, P, &Y[0] );
	  
	  cp = _cantera_gas.cp_mass();
	}
      catch(Cantera::CanteraError)
	{
	  Cantera::showErrors(std::cerr);
	  libmesh_error();
	}

    }

    return cp;
  }

  libMesh::Real CanteraThermodynamics::cv( const CachedValues& cache, unsigned int qp ) const
  {
    libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const libMesh::Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];
    
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    libMesh::Real cv = 0.0;

    {
      //libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

      /*! \todo Need to make sure this will work in a threaded environment.
	Not sure if we will get thread lock here or not. */
      try
	{
	  _cantera_gas.setState_TPY( T, P, &Y[0] );
	  
	  cv = _cantera_gas.cv_mass();
	}
      catch(Cantera::CanteraError)
	{
	  Cantera::showErrors(std::cerr);
	  libmesh_error();
	}

    }

    return cv;
  }

  libMesh::Real CanteraThermodynamics::h( const CachedValues& cache, unsigned int qp,
				 unsigned int species ) const
  {
    libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const libMesh::Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    std::vector<libMesh::Real> h_RT( Y.size(), 0.0 );

    {
      //libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);
    
      /*! \todo Need to make sure this will work in a threaded environment.
	Not sure if we will get thread lock here or not. */
      try
	{
	  _cantera_gas.setState_TPY( T, P, &Y[0] );
	  
	  _cantera_gas.getEnthalpy_RT( &h_RT[0] );
	}
      catch(Cantera::CanteraError)
	{
	  Cantera::showErrors(std::cerr);
	  libmesh_error();
	}

    }

    return h_RT[species]*_chem_mixture.R(species)*T;
  }

  void CanteraThermodynamics::h( const CachedValues& cache, unsigned int qp,
				 std::vector<libMesh::Real>& h) const
  {
    libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const libMesh::Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    libmesh_assert_equal_to( Y.size(), h.size() );
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    {
      //libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);
    
      /*! \todo Need to make sure this will work in a threaded environment.
	Not sure if we will get thread lock here or not. */
      try
	{
	  _cantera_gas.setState_TPY( T, P, &Y[0] );
	  
	  _cantera_gas.getEnthalpy_RT( &h[0] );
	}
      catch(Cantera::CanteraError)
	{
	  Cantera::showErrors(std::cerr);
	  libmesh_error();
	}

    for( unsigned int s = 0; s < h.size(); s++ )
      {
	h[s] *= _chem_mixture.R(s)*T;
      }

    }

    return;
  }

} // namespace GRINS

#endif //GRINS_HAVE_CANTERA
