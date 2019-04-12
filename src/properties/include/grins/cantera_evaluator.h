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
  //! Wrapper class for evaluating thermochemistry and transport properties using Cantera
  /*!
    This class is expected to be constructed *after* threads have been forked and will only
    live during the lifetime of the thread. Note that this documentation will always
    be built regardless if Cantera is included in the GRINS build or not. Check configure
    output to confirm that Cantera was included in the build if you wish to use it.
  */
  class CanteraEvaluator
  {
  public:

    CanteraEvaluator( CanteraMixture& mixture );
    ~CanteraEvaluator(){};

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
    libMesh::Real cp( const libMesh::Real& T, const libMesh::Real P, const std::vector<libMesh::Real>& Y );
    
    void cp_s( const libMesh::Real& T, const libMesh::Real P, const std::vector<libMesh::Real>& Y, std::vector<libMesh::Real>& Cp_s);

    libMesh::Real cv( const libMesh::Real& T, const libMesh::Real P, const std::vector<libMesh::Real>& Y );

    libMesh::Real h_s( const libMesh::Real& T, unsigned int species );

    // Transport
    libMesh::Real mu( const libMesh::Real& T, const libMesh::Real P, const std::vector<libMesh::Real>& Y );

    libMesh::Real k( const libMesh::Real& T, const libMesh::Real P, const std::vector<libMesh::Real>& Y );

    void mu_and_k( const libMesh::Real& T, const libMesh::Real P, const std::vector<libMesh::Real>& Y,
                   libMesh::Real& mu, libMesh::Real& k );

    void mu_and_k_and_D( const libMesh::Real T,
                         const libMesh::Real rho,
                         const libMesh::Real cp,
                         const std::vector<libMesh::Real>& Y,
                         libMesh::Real& mu, libMesh::Real& k,
                         std::vector<libMesh::Real>& D );

    // Kinetics
    void omega_dot( const libMesh::Real& T, libMesh::Real rho,
                    const std::vector<libMesh::Real> mass_fractions,
                    std::vector<libMesh::Real>& omega_dot );

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
  libMesh::Real CanteraEvaluator::cp( const libMesh::Real& T,
                                      const libMesh::Real P,
                                      const std::vector<libMesh::Real>& Y )
  {
    return _thermo.cp(T,P,Y);
  }

  inline
  void CanteraEvaluator::cp_s( const libMesh::Real& T,
				 const libMesh::Real P,
				 const std::vector<libMesh::Real>& Y,
				 std::vector<libMesh::Real>& Cp_s)
  {
    _thermo.cp_s(T,P,Y,Cp_s);
    return;
  }

  inline
  libMesh::Real CanteraEvaluator::cv( const libMesh::Real& T,
                                      const libMesh::Real P,
                                      const std::vector<libMesh::Real>& Y )
  {
    return _thermo.cv(T,P,Y);
  }

  inline
  libMesh::Real CanteraEvaluator::h_s( const libMesh::Real& T, unsigned int species )
  {
    return _thermo.h(T,species);
  }

  inline
  libMesh::Real CanteraEvaluator::mu( const libMesh::Real& T,
                                      const libMesh::Real P,
                                      const std::vector<libMesh::Real>& Y )
  {
    return _transport.mu(T,P,Y);
  }

  inline
  libMesh::Real CanteraEvaluator::k( const libMesh::Real& T,
                                     const libMesh::Real P,
                                     const std::vector<libMesh::Real>& Y )
  {
    return _transport.k(T,P,Y);
  }

  inline
  void CanteraEvaluator::mu_and_k( const libMesh::Real& T,
                                   const libMesh::Real P,
                                   const std::vector<libMesh::Real>& Y,
                                   libMesh::Real& mu, libMesh::Real& k )
  {
    mu = _transport.mu(T,P,Y);
    k = _transport.k(T,P,Y);
  }

  inline
  void CanteraEvaluator::mu_and_k_and_D( const libMesh::Real T,
                                         const libMesh::Real rho,
                                         const libMesh::Real cp,
                                         const std::vector<libMesh::Real>& Y,
                                         libMesh::Real& mu, libMesh::Real& k,
                                         std::vector<libMesh::Real>& D )
  {
    _transport.mu_and_k_and_D( T, rho, cp, Y, mu, k, D );
  }

  inline
  void CanteraEvaluator::omega_dot( const libMesh::Real& T, libMesh::Real rho,
                                    const std::vector<libMesh::Real> mass_fractions,
                                    std::vector<libMesh::Real>& omega_dot )
  {
    _kinetics.omega_dot(T,rho,mass_fractions,omega_dot);
  }

} // end namespace GRINS

#endif // GRINS_HAVE_CANTERA

#endif // GRINS_CANTERA_EVALUATOR_H
