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

#ifndef GRINS_MIXTURE_TRANSPORT_H
#define GRINS_MIXTURE_TRANSPORT_H

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/libmesh_common.h"

// GRINS
#include "grins/chemical_mixture.h"

namespace GRINS
{
  template<typename Viscosity, typename Conductivity, typename Diffusivity>
  class MixtureTransport
  {
  public:

    MixtureTransport( const GetPot& input, const ChemicalMixture& chem_mixture );
    ~MixtureTransport();

    libMesh::Real mu( libMesh::Real T, unsigned int species );
    libMesh::Real mu( libMesh::Real T );

  protected:
    
    const ChemicalMixture& _chem_mixture;

    Viscosity _viscosity;

    Conductivity _thermal_conductivity;

    Diffusivity _diffusivity;

  };

} // namespace GRINS

#endif //GRINS_MIXTURE_TRANSPORT_H
