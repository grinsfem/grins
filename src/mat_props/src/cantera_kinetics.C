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

#include "cantera_kinetics.h"

#ifdef HAVE_CANTERA

namespace GRINS
{
  CanteraKinetics::CanteraKinetics( const GetPot& input, const ChemicalMixture& chem_mixture )
    : _chem_mixture(chem_mixture),
      _cantera_gas( CanteraSingleton::cantera_instance(input) )
  {
    return;
  }

  CanteraKinetics::~CanteraKinetics()
  {
    return;
  }

  void CanteraKinetics::omega_dot( const ReactingFlowCache& cache, std::vector<Real>& omega_dot )
  {
    const Real T = cache.T();
    const Real P = cache.P();
    const std::vector<Real>& Y = cache.mass_fractions();

    libmesh_assert_equal_to( Y.size(), omega_dot.size() );

    Threads::spin_mutex cantera_mutex;
    Threads::spin_mutex::scoped_lock lock(cantera_mutex);
    
    /*! \todo Need to make sure this will work in a threaded environment.
      Not sure if we will get thread lock here or not. */
    try
      {
	_cantera_gas.setState_TPY(T, P, &Y[0]);
	_cantera_gas.getNetProductionRates(&omega_dot[0]);
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
