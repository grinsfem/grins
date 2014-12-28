//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2014 Paul T. Bauman, Roy H. Stogner
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


#ifndef GRINS_CANTERA_EVALUATOR_H
#define GRINS_CANTERA_EVALUATOR_H

#include "grins_config.h"

#ifdef GRINS_HAVE_CANTERA

// GRINS
#include "grins/cantera_mixture.h"
#include "grins/cantera_thermo.h"
#include "grins/cantera_transport.h"
#include "grins/cantera_kinetics.h"

namespace GRINS
{
  class CanteraEvaluator
  {
  public:

    CanteraEvaluator( CanteraMixture& mixture );
    ~CanteraEvaluator();

    // Chemistry
    libMesh::Real M( unsigned int species ) const;

    libMesh::Real M_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real R( unsigned int species ) const;

    libMesh::Real R_mix( const std::vector<libMesh::Real>& mass_fractions ) const;

    libMesh::Real X( unsigned int species, libMesh::Real M, libMesh::Real mass_fraction ) const;

    void X( libMesh::Real M, const std::vector<libMesh::Real>& mass_fractions, 
	    std::vector<libMesh::Real>& mole_fractions ) const;

    unsigned int species_index( const std::string& species_name ) const;

    std::string species_name( unsigned int species_index ) const;

    // Thermo
    libMesh::Real cp( const CachedValues& cache, unsigned int qp ) const;

    libMesh::Real cv( const CachedValues& cache, unsigned int qp ) const;

    libMesh::Real h_s(const CachedValues& cache, unsigned int qp, unsigned int species) const;

    void h_s(const CachedValues& cache, unsigned int qp, std::vector<libMesh::Real>& h) const;

    libMesh::Real h_s( const libMesh::Real& T, unsigned int species );

    // Transport
    libMesh::Real mu( const CachedValues& cache, unsigned int qp ) const;

    libMesh::Real k( const CachedValues& cache, unsigned int qp ) const;

    void mu_and_k( const CachedValues& cache, unsigned int qp,
                   libMesh::Real& mu, libMesh::Real& k );

    void D( const CachedValues& cache, unsigned int qp,
	    std::vector<libMesh::Real>& D ) const;

    // Kinetics
    void omega_dot( const CachedValues& cache, unsigned int qp,
		    std::vector<libMesh::Real>& omega_dot ) const;

    libMesh::Real cp( const libMesh::Real& /*T*/,
                      const std::vector<libMesh::Real>& /*Y*/ )
    {
      libmesh_not_implemented();
      return 0.0;
    }

    libMesh::Real mu( const libMesh::Real& /*T*/,
                      const std::vector<libMesh::Real>& /*Y*/ )
    {
      libmesh_not_implemented();
      return 0.0;
    }

    libMesh::Real k( const libMesh::Real& /*T*/,
                     const std::vector<libMesh::Real>& /*Y*/ )
    {
      libmesh_not_implemented();
      return 0.0;
    }

    void D( const libMesh::Real /*rho*/,
            const libMesh::Real /*cp*/,
            const libMesh::Real /*k*/,
	    std::vector<libMesh::Real>& /*D*/ )
    {
      libmesh_not_implemented();
      return;
    }

  protected:

    CanteraMixture& _chem;

    CanteraThermodynamics _thermo;

    CanteraTransport _transport;

    CanteraKinetics _kinetics;

  private:

    CanteraEvaluator();

  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  libMesh::Real CanteraEvaluator::M( unsigned int species ) const
  {
    return _chem.M(species);
  }

  inline
  libMesh::Real CanteraEvaluator::M_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    return _chem.M_mix(mass_fractions);
  }

  inline
  libMesh::Real CanteraEvaluator::R( unsigned int species ) const
  {
    return _chem.R(species);
  }

  inline
  libMesh::Real CanteraEvaluator::R_mix( const std::vector<libMesh::Real>& mass_fractions ) const
  {
    return _chem.R_mix(mass_fractions);
  }
  
  inline
  libMesh::Real CanteraEvaluator::X( unsigned int species, libMesh::Real M, libMesh::Real mass_fraction ) const
  {
    return _chem.X(species,M,mass_fraction);
  }
  
  inline
  void CanteraEvaluator::X( libMesh::Real M, const std::vector<libMesh::Real>& mass_fractions, 
                            std::vector<libMesh::Real>& mole_fractions ) const
  {
    _chem.X(M,mass_fractions,mole_fractions);
    return;
  }
  
  inline
  unsigned int CanteraEvaluator::species_index( const std::string& species_name ) const
  {
    return _chem.species_index(species_name);
  }
  
  inline
  std::string CanteraEvaluator::species_name( unsigned int species_index ) const
  {
    return _chem.species_name(species_index);
  }

  inline
  libMesh::Real CanteraEvaluator::cp( const CachedValues& cache, unsigned int qp ) const
  {
    return _thermo.cp(cache,qp);
  }

  inline
  libMesh::Real CanteraEvaluator::cv( const CachedValues& cache, unsigned int qp ) const
  {
    return _thermo.cv(cache,qp);
  }
  
  inline
  libMesh::Real CanteraEvaluator::h_s(const CachedValues& cache, unsigned int qp, unsigned int species) const
  {
    return _thermo.h(cache,qp,species);
  }
  
  inline
  void CanteraEvaluator::h_s(const CachedValues& cache, unsigned int qp, std::vector<libMesh::Real>& h) const
  {
    _thermo.h(cache,qp,h);
    return;
  }

  inline
  libMesh::Real CanteraEvaluator::h_s( const libMesh::Real& T, unsigned int species )
  {
    return _thermo.h(T,species);
  }

  inline
  libMesh::Real CanteraEvaluator::mu( const CachedValues& cache, unsigned int qp ) const
  {
    return _transport.mu(cache,qp);
  }

  inline
  libMesh::Real CanteraEvaluator::k( const CachedValues& cache, unsigned int qp ) const
  {
    return _transport.k(cache,qp);
  }

  inline
  void CanteraEvaluator::mu_and_k( const CachedValues& cache, unsigned int qp,
                                   libMesh::Real& mu, libMesh::Real& k )
  {
    mu = _transport.mu(cache,qp);
    k = _transport.k(cache,qp);
    return;
  }

  inline
  void CanteraEvaluator::D( const CachedValues& cache, unsigned int qp,
                            std::vector<libMesh::Real>& D ) const
  {
    return _transport.D(cache,qp,D);
  }

  inline
  void CanteraEvaluator::omega_dot( const CachedValues& cache, unsigned int qp,
                                    std::vector<libMesh::Real>& omega_dot ) const
  {
    return _kinetics.omega_dot(cache,qp,omega_dot);
  }

} // end namespace GRINS

#endif // GRINS_HAVE_CANTERA

#endif // GRINS_CANTERA_EVALUATOR_H
