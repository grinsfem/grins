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


#ifndef GRINS_ANTIOCH_MIXTURE_H
#define GRINS_ANTIOCH_MIXTURE_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// GRINS
#include "grins/antioch_chemistry.h"
#include "grins/property_types.h"

// libMesh
#include "libmesh/libmesh_common.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/chemical_mixture.h"
#include "antioch/cea_mixture.h"
#include "antioch/reaction_set.h"

// Boost
#include "boost/scoped_ptr.hpp"

// libMesh forward declarations
class GetPot;

namespace GRINS
{
  class AntiochMixture : public AntiochChemistry
  {
  public:

    AntiochMixture( const GetPot& input );

    virtual ~AntiochMixture();

    const Antioch::ReactionSet<libMesh::Real>& reaction_set() const;

    const Antioch::CEAThermoMixture<libMesh::Real>& cea_mixture() const;

    libMesh::Real h_stat_mech_ref_correction( unsigned int species ) const;

  protected:

    boost::scoped_ptr<Antioch::ReactionSet<libMesh::Real> > _reaction_set;

    boost::scoped_ptr<Antioch::CEAThermoMixture<libMesh::Real> > _cea_mixture;

    std::vector<libMesh::Real> _h_stat_mech_ref_correction;

    void build_stat_mech_ref_correction();

  private:

    AntiochMixture();

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  const Antioch::ReactionSet<libMesh::Real>& AntiochMixture::reaction_set() const
  {
    return *_reaction_set.get();
  }

  inline
  const Antioch::CEAThermoMixture<libMesh::Real>& AntiochMixture::cea_mixture() const
  {
    return *_cea_mixture.get();
  }

  inline
  libMesh::Real AntiochMixture::h_stat_mech_ref_correction( unsigned int species ) const
  {
    return _h_stat_mech_ref_correction[species];
  }
  
} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_MIXTURE_H
