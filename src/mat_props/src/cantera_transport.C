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

#include "cantera_transport.h"

#ifdef HAVE_CANTERA

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

  Real CanteraTransport::mu( const ReactingFlowCache& cache )
  {
    const Real T = cache.T();
    const Real P = cache.P();
    const std::vector<Real>& Y = cache.mass_fractions();

    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    Real mu = 0.0;

    Threads::spin_mutex cantera_mutex;
    Threads::spin_mutex::scoped_lock lock(cantera_mutex);
    
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

    lock.release();

    return mu;
  }

  Real CanteraTransport::k( const ReactingFlowCache& cache )
  {
    const Real T = cache.T();
    const Real P = cache.P();
    const std::vector<Real>& Y = cache.mass_fractions();

    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    Real k = 0.0;

    Threads::spin_mutex cantera_mutex;
    Threads::spin_mutex::scoped_lock lock(cantera_mutex);
    
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

    lock.release();

    return k;
  }

  void CanteraTransport::D( const ReactingFlowCache& cache, std::vector<Real>& D )
  {
    const Real T = cache.T();
    const Real P = cache.P();
    const std::vector<Real>& Y = cache.mass_fractions();

    libmesh_assert_equal_to( Y.size(), D.size() );
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    Threads::spin_mutex cantera_mutex;
    Threads::spin_mutex::scoped_lock lock(cantera_mutex);
    
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

    lock.release();

    return;
  }

} // namespace GRINS

#endif //HAVE_CANTERA
