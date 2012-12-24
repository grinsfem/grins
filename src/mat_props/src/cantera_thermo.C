//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "cantera_thermo.h"

#ifdef GRINS_HAVE_CANTERA

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

  Real CanteraThermodynamics::cp( const ReactingFlowCache& cache )
  {
    const Real T = cache.T();
    
    const Real P = cache.P();

    const std::vector<Real>& Y = cache.mass_fractions();
    
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    Real cp = 0.0;
 
    Threads::spin_mutex cantera_mutex;
    Threads::spin_mutex::scoped_lock lock(cantera_mutex);

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

    lock.release();

    return cp;
  }

  Real CanteraThermodynamics::cp( const CachedValues& cache, unsigned int qp )
  {
    const Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];
    
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    Real cp = 0.0;

    {
      Threads::spin_mutex cantera_mutex;
      Threads::spin_mutex::scoped_lock lock(cantera_mutex);

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

      lock.release();
    }

    return cp;
  }

  Real CanteraThermodynamics::cv( const ReactingFlowCache& cache )
  {
    const Real T = cache.T();
    
    const Real P = cache.P();

    const std::vector<Real>& Y = cache.mass_fractions();
    
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    Real cv = 0.0;

    Threads::spin_mutex cantera_mutex;
    Threads::spin_mutex::scoped_lock lock(cantera_mutex);

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

    lock.release();

    return cv;
  }

  Real CanteraThermodynamics::cv( const CachedValues& cache, unsigned int qp )
  {
    const Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];
    
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    Real cv = 0.0;

    Threads::spin_mutex cantera_mutex;
    Threads::spin_mutex::scoped_lock lock(cantera_mutex);

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

    lock.release();

    return cv;
  }

  Real CanteraThermodynamics::h(const ReactingFlowCache& cache, unsigned int species)
  {
    const Real T = cache.T();
    
    const Real P = cache.P();

    const std::vector<Real>& Y = cache.mass_fractions();

    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    std::vector<Real> h_RT( Y.size(), 0.0 );

    Threads::spin_mutex cantera_mutex;
    Threads::spin_mutex::scoped_lock lock(cantera_mutex);
    
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

    lock.release();

    return h_RT[species]*_chem_mixture.R(species)*T;
  }

  Real CanteraThermodynamics::h(const CachedValues& cache, unsigned int qp,
				unsigned int species)
  {
    const Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    std::vector<Real> h_RT( Y.size(), 0.0 );

    {
      Threads::spin_mutex cantera_mutex;
      Threads::spin_mutex::scoped_lock lock(cantera_mutex);
    
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

      lock.release();
    }

    return h_RT[species]*_chem_mixture.R(species)*T;
  }

  void CanteraThermodynamics::h(const ReactingFlowCache& cache, std::vector<Real>& h)
  {
    const Real T = cache.T();
    
    const Real P = cache.P();

    const std::vector<Real>& Y = cache.mass_fractions();

    libmesh_assert_equal_to( Y.size(), h.size() );
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    Threads::spin_mutex cantera_mutex;
    Threads::spin_mutex::scoped_lock lock(cantera_mutex);
    
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

    lock.release();

    for( unsigned int s = 0; s < h.size(); s++ )
      {
	h[s] *= _chem_mixture.R(s)*T;
      }

    return;
  }

  void CanteraThermodynamics::h( const CachedValues& cache, unsigned int qp,
				 std::vector<Real>& h)
  {
    const Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    libmesh_assert_equal_to( Y.size(), h.size() );
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    {
      Threads::spin_mutex cantera_mutex;
      Threads::spin_mutex::scoped_lock lock(cantera_mutex);
    
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

      lock.release();
    }

    for( unsigned int s = 0; s < h.size(); s++ )
      {
	h[s] *= _chem_mixture.R(s)*T;
      }

    return;
  }

} // namespace GRINS

#endif //GRINS_HAVE_CANTERA
