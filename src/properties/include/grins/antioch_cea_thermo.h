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

#ifndef GRINS_ANTIOCH_CEA_THERMO_H
#define GRINS_ANTIOCH_CEA_THERMO_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// libMesh
#include "libmesh/libmesh_common.h"

// Antioch
#include "antioch/vector_utils_decl.h"
#include "antioch/vector_utils.h"
#include "antioch/temp_cache.h"
#include "antioch/cea_evaluator.h"

namespace GRINS
{
  // GRINS forward declarations
  class AntiochMixture;

  class AntiochCEAThermo
  {
  public:

    AntiochCEAThermo( const AntiochMixture& mixture );
    ~AntiochCEAThermo();

    libMesh::Real cp( const Antioch::TempCache<libMesh::Real>& cache, unsigned int species ) const;

    libMesh::Real cp( const Antioch::TempCache<libMesh::Real>& cache,
                      const std::vector<libMesh::Real> mass_fractions ) const;

    libMesh::Real cv( const Antioch::TempCache<libMesh::Real>& cache, unsigned int species ) const;

    libMesh::Real cv( const Antioch::TempCache<libMesh::Real>& cache,
                      const std::vector<libMesh::Real> mass_fractions ) const;

    libMesh::Real h_s( const Antioch::TempCache<libMesh::Real>& cache, unsigned int species ) const;

    void h_s( const Antioch::TempCache<libMesh::Real>& cache, std::vector<libMesh::Real>& h_s ) const;

  protected:
    
    const AntiochMixture& _antioch_mixture;

    const Antioch::CEAEvaluator<libMesh::Real>& _antioch_cea_evaluator;

  private:

    AntiochCEAThermo();
    
  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  libMesh::Real AntiochCEAThermo::cp( const Antioch::TempCache<libMesh::Real>& cache, unsigned int species ) const
  {
    return _antioch_cea_evaluator.cp( cache, species );
  }

  inline
  libMesh::Real AntiochCEAThermo::cp( const Antioch::TempCache<libMesh::Real>& cache,
                                      const std::vector<libMesh::Real> mass_fractions ) const
  {
    return _antioch_cea_evaluator.cp( cache, mass_fractions );
  }

  inline
  libMesh::Real AntiochCEAThermo::cv( const Antioch::TempCache<libMesh::Real>& cache, unsigned int species ) const
  {
    return _antioch_cea_evaluator.cv( cache, species );
  }

  inline
  libMesh::Real AntiochCEAThermo::cv( const Antioch::TempCache<libMesh::Real>& cache,
                                      const std::vector<libMesh::Real> mass_fractions ) const
  {
    return _antioch_cea_evaluator.cv( cache, mass_fractions );
  }

  inline
  libMesh::Real AntiochCEAThermo::h_s( const Antioch::TempCache<libMesh::Real>& cache, unsigned int species ) const
  {
    return _antioch_cea_evaluator.h( cache, species );
  }

  inline
  void AntiochCEAThermo::h_s( const Antioch::TempCache<libMesh::Real>& cache,
                              std::vector<libMesh::Real>& h_s ) const
  {
    _antioch_cea_evaluator.h( cache, h_s );
    return;
  }

} // end namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_ANTIOCH_CEA_THERMO_H
