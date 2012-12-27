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

#ifndef GRINS_ARRHENIUS_RATE_H
#define GRINS_ARRHENIUS_RATE_H

// libMesh
#include "libmesh_common.h"

namespace GRINS
{
  //! Modified Arrhenius rate equation.
  /*!
   * Modified Arrhenius rate equation.  Computes rates of the form
   * \f$ C_f\times T^\eta\times \exp(-E_a/T) \f$. This class copied from
   * the \p FIN-S code and slightly reformatted for \p GRINS.
   */
  class ArrheniusRate
  {
  
  public:

    ArrheniusRate (const Real Cf=0., const Real eta=0., const Real Ea=0.);
    ~ArrheniusRate();
    
    void set_Cf( const Real Cf );
    void set_eta( const Real eta );
    void set_Ea( const Real Ea );

    Real Cf() const;
    Real eta() const;
    Real Ea() const;

    //! \return the rate evaluated at \p T.
    Real operator()(const Real T) const;

    //! \return the derivative with respect to temperature evaluated at \p T.
    Real derivative( const Real T ) const;

    //! Simultaneously evaluate the rate and its derivative at \p T.
    void rate_and_derivative(const Real T, Real& rate, Real& drate_dT) const;

    //! Formatted print, by default to \p libMesh::out.
    void print(std::ostream& os = libMesh::out) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const ArrheniusRate& rate);

  private:

    Real _Cf;
    Real _eta;
    Real _Ea;
    
  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  void ArrheniusRate::set_Cf( const Real Cf )
  {
    _Cf = Cf;
    return;
  }

  inline
  void ArrheniusRate::set_eta( const Real eta )
  {
    _eta = eta;
    return;
  }

  inline
  void ArrheniusRate::set_Ea( const Real Ea )
  {
    _Ea = Ea;
    return;
  }

  inline
  Real ArrheniusRate::Cf() const
  { return _Cf; }

  inline
  Real ArrheniusRate::eta() const
  { return _eta; }

  inline
  Real ArrheniusRate::Ea() const
  { return _Ea; }

  inline
  Real ArrheniusRate::operator()(const Real T) const
  {
    return _Cf* (std::pow(T,_eta)*std::exp(-_Ea/T));
  }

  inline
  Real ArrheniusRate::derivative( const Real T ) const
  {
    return (*this)(T)/T*(_eta + _Ea/T);
  }

  inline
  void ArrheniusRate::rate_and_derivative (const Real T,
					   Real& rate,
					   Real& drate_dT) const
    {
      rate     = (*this)(T);
      drate_dT = rate/T*(_eta + _Ea/T);
      return;
    }

} // namespace GRINS

#endif // GRINS_ARRHENIUS_RATE_H
