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
#include "libmesh/libmesh_common.h"

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

    ArrheniusRate (const libMesh::Real Cf=0., const libMesh::Real eta=0., const libMesh::Real Ea=0.);
    ~ArrheniusRate();
    
    void set_Cf( const libMesh::Real Cf );
    void set_eta( const libMesh::Real eta );
    void set_Ea( const libMesh::Real Ea );

    void scale_Ea( const libMesh::Real scale );

    libMesh::Real Cf() const;
    libMesh::Real eta() const;
    libMesh::Real Ea() const;

    //! \return the rate evaluated at \p T.
    libMesh::Real operator()(const libMesh::Real T) const;

    //! \return the derivative with respect to temperature evaluated at \p T.
    libMesh::Real derivative( const libMesh::Real T ) const;

    //! Simultaneously evaluate the rate and its derivative at \p T.
    void rate_and_derivative(const libMesh::Real T, libMesh::Real& rate, libMesh::Real& drate_dT) const;

    //! Formatted print, by default to \p libMesh::out.
    void print(std::ostream& os = libMesh::out) const;

    //! Formatted print.
    friend std::ostream& operator<<(std::ostream& os, const ArrheniusRate& rate);

  private:

    libMesh::Real _Cf;
    libMesh::Real _eta;
    libMesh::Real _Ea;
    
  };

  /* ------------------------- Inline Functions -------------------------*/
  inline
  void ArrheniusRate::set_Cf( const libMesh::Real Cf )
  {
    _Cf = Cf;
    return;
  }

  inline
  void ArrheniusRate::set_eta( const libMesh::Real eta )
  {
    _eta = eta;
    return;
  }

  inline
  void ArrheniusRate::set_Ea( const libMesh::Real Ea )
  {
    _Ea = Ea;
    return;
  }

  inline
  void ArrheniusRate::scale_Ea( const libMesh::Real scale )
  {
    _Ea *= scale;
    return;
  }

  inline
  libMesh::Real ArrheniusRate::Cf() const
  { return _Cf; }

  inline
  libMesh::Real ArrheniusRate::eta() const
  { return _eta; }

  inline
  libMesh::Real ArrheniusRate::Ea() const
  { return _Ea; }

  inline
  libMesh::Real ArrheniusRate::operator()(const libMesh::Real T) const
  {
    return _Cf* (std::pow(T,_eta)*std::exp(-_Ea/T));
  }

  inline
  libMesh::Real ArrheniusRate::derivative( const libMesh::Real T ) const
  {
    return (*this)(T)/T*(_eta + _Ea/T);
  }

  inline
  void ArrheniusRate::rate_and_derivative (const libMesh::Real T,
					   libMesh::Real& rate,
					   libMesh::Real& drate_dT) const
    {
      rate     = (*this)(T);
      drate_dT = rate/T*(_eta + _Ea/T);
      return;
    }

} // namespace GRINS

#endif // GRINS_ARRHENIUS_RATE_H
