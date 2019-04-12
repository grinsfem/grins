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

#ifdef GRINS_HAVE_CANTERA

// This class
#include "grins/cantera_thermo.h"

// GRINS
#include "grins/cantera_mixture.h"
#include "grins/cached_values.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{

  CanteraThermodynamics::CanteraThermodynamics( CanteraMixture& mixture )
    : _cantera_mixture(mixture),
      _cantera_gas(mixture.get_chemistry())
  {}

  libMesh::Real CanteraThermodynamics::cp( const libMesh::Real& T,
                                           const libMesh::Real P,
                                           const std::vector<libMesh::Real>& Y )
  {
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    libMesh::Real cp = 0.0;

    {
      /*! \todo Need to make sure this will work in a threaded environment.
        Not sure if we will get thread lock here or not. */
      libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

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

  void CanteraThermodynamics::cp_s( const libMesh::Real& T,
				    const libMesh::Real P,
				    const std::vector<libMesh::Real>& Y,
				    std::vector<libMesh::Real> &Cp_s)
  {
    libmesh_assert_equal_to( Cp_s.size(), _cantera_gas.nSpecies() );
    {
      libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);
      
      try
	{
	  _cantera_gas.setState_TPY( T, P, &Y[0]);
	  _cantera_gas.getPartialMolarCp(&Cp_s[0]);
	  for(unsigned int species = 0; species < _cantera_gas.nSpecies(); species++)
	    Cp_s[species] = Cp_s[species]/_cantera_gas.molecularWeight(species);
	}
      catch(Cantera::CanteraError)
	{
	  Cantera::showErrors(std::cerr);
	  libmesh_error();
	}
    }
    
    return;
  }

  libMesh::Real CanteraThermodynamics::cv( const libMesh::Real& T,
                                           const libMesh::Real P,
                                           const std::vector<libMesh::Real>& Y )
  {
    libmesh_assert_equal_to( Y.size(), _cantera_gas.nSpecies() );

    libMesh::Real cv = 0.0;

    {
      /*! \todo Need to make sure this will work in a threaded environment.
        Not sure if we will get thread lock here or not. */
      libMesh::Threads::spin_mutex::scoped_lock lock(cantera_mutex);

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

  libMesh::Real CanteraThermodynamics::h( const libMesh::Real& T, unsigned int species )
  {
    std::vector<libMesh::Real> h_RT( _cantera_gas.nSpecies(), 0.0 );

    try
      {
        _cantera_gas.setTemperature( T );

        _cantera_gas.getEnthalpy_RT( &h_RT[0] );
      }
    catch(Cantera::CanteraError)
      {
        Cantera::showErrors(std::cerr);
        libmesh_error();
      }

    return h_RT[species]*_cantera_mixture.R(species)*T;
  }

} // namespace GRINS

#endif //GRINS_HAVE_CANTERA
