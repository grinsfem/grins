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

#ifdef GRINS_HAVE_CANTERA

// This class
#include "grins/cantera_mixture.h"

// libMesh
#include "libmesh/getpot.h"

namespace GRINS
{
  CanteraMixture::CanteraMixture( const GetPot& input )
    : _cantera_gas(NULL),
      _cantera_transport(NULL)
  {
    const std::string cantera_chem_file = input( "Physics/Chemistry/chem_file", "DIE!" );
    const std::string mixture = input( "Physics/Chemistry/mixture", "DIE!" );

    try
      {
        _cantera_gas.reset( new Cantera::IdealGasMix( cantera_chem_file, mixture ) );
      }
    catch(Cantera::CanteraError)
      {
        Cantera::showErrors(std::cerr);
        libmesh_error();
      }

    try
      {
        _cantera_transport.reset( Cantera::newTransportMgr("Mix", _cantera_gas.get()) );
      }
    catch(Cantera::CanteraError)
      {
        Cantera::showErrors(std::cerr);
        libmesh_error();
      }

    return;
  }

  CanteraMixture::~CanteraMixture()
  {
    return;
  }

  libMesh::Real CanteraMixture::M_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    libmesh_assert_equal_to( mass_fractions.size(), _cantera_gas->nSpecies() );
    
    libMesh::Real M = 0;
    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      {
	M += mass_fractions[s]/(this->M(s));
      }

    return 1.0/M;
  }

  libMesh::Real CanteraMixture::R_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    libmesh_assert_equal_to( mass_fractions.size(), _cantera_gas->nSpecies() );
    
    libMesh::Real R = 0.0;
    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      {
	R += mass_fractions[s]*this->R(s);
      }
    
    return R;
  }

  void CanteraMixture::X( libMesh::Real M_mix, const std::vector<libMesh::Real>& mass_fractions, 
                          std::vector<libMesh::Real>& mole_fractions ) const
  {
    libmesh_assert_equal_to( mass_fractions.size(), _cantera_gas->nSpecies() );

    libmesh_assert_equal_to( mole_fractions.size(), mass_fractions.size() );

    for( unsigned int s = 0; s < mass_fractions.size(); s++ )
      {
	mole_fractions[s] = this->X(s, M_mix, mass_fractions[s]);
      }

    return;
  }

}// end namespace GRINS

#endif // GRINS_HAVE_CANTERA
