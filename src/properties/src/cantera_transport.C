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

#include "grins/cantera_transport.h"

#ifdef GRINS_HAVE_CANTERA

namespace
{
  libMesh::Threads::spin_mutex transport_mutex;
}

namespace GRINS
{

  CanteraTransport::CanteraTransport( const GetPot& input, const ChemicalMixture& chem_mixture )
    : _chem_mixture(chem_mixture),
      _cantera_gas( CanteraSingleton::cantera_instance(input) ),
      _cantera_transport( Cantera::newTransportMgr("Mix", &_cantera_gas) )
  {
    return;
  }

  CanteraTransport::~CanteraTransport()
  {
    delete _cantera_transport;
    return;
  }

  libMesh::Real CanteraTransport::mu( const CachedValues& cache, unsigned int qp ) const
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const libMesh::Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    libMesh::Real mu = 0.0;

    {
      libMesh::Threads::spin_mutex::scoped_lock lock(transport_mutex);
    
      /*! \todo Need to make sure this will work in a threaded environment.
	Not sure if we will get thread lock here or not. */
      try
	{
	  _cantera_gas.setState_TPY(T, P, &Y[0]);
	  mu =  _cantera_transport->viscosity();
	}
      catch(Cantera::CanteraError)
	{
	  Cantera::showErrors(std::cerr);
	  libmesh_error();
	}

    }

    return mu;
  }

  libMesh::Real CanteraTransport::k( const CachedValues& cache, unsigned int qp ) const
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const libMesh::Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    libMesh::Real k = 0.0;

    {
      libMesh::Threads::spin_mutex::scoped_lock lock(transport_mutex);
    
      /*! \todo Need to make sure this will work in a threaded environment.
	Not sure if we will get thread lock here or not. */
      try
	{
	  _cantera_gas.setState_TPY(T, P, &Y[0]);
	  k =  _cantera_transport->thermalConductivity();
	}
      catch(Cantera::CanteraError)
	{
	  Cantera::showErrors(std::cerr);
	  libmesh_error();
	}

    }

    return k;
  }

  void CanteraTransport::D( const CachedValues& cache, unsigned int qp,
			    std::vector<libMesh::Real>& D ) const
  {
    const libMesh::Real T = cache.get_cached_values(Cache::TEMPERATURE)[qp];
    const libMesh::Real P = cache.get_cached_values(Cache::THERMO_PRESSURE)[qp];
    const std::vector<libMesh::Real>& Y = cache.get_cached_vector_values(Cache::MASS_FRACTIONS)[qp];

    libmesh_assert_equal_to( Y.size(), D.size() );
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    {
      libMesh::Threads::spin_mutex::scoped_lock lock(transport_mutex);
    
      /*! \todo Need to make sure this will work in a threaded environment.
	Not sure if we will get thread lock here or not. */
      try
	{
	  _cantera_gas.setState_TPY(T, P, &Y[0]);
	  _cantera_transport->getMixDiffCoeffsMass(&D[0]);
	}
      catch(Cantera::CanteraError)
	{
	  Cantera::showErrors(std::cerr);
	  libmesh_error();
	}

    }

    return;
  }

} // namespace GRINS

#endif //GRINS_HAVE_CANTERA
