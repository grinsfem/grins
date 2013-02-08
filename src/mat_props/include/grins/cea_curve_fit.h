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

#ifndef GRINS_CEA_CURVE_FIT_H
#define GRINS_CEA_CURVE_FIT_H

// GRINS
#include "grins/chemical_mixture.h"

namespace GRINS
{
  class CEACurveFit
  {
  public:
    
    CEACurveFit( const std::vector<libMesh::Real>& coeffs );
    ~CEACurveFit();

    //! The number of intervals for this CEA curve fit
    inline
    unsigned int n_intervals() const 
    { return _coefficients.size() / _n_coeffs; }

    //! The interval the input temperature lies in
    /*!
      @returns which curve fit interval the input temperature 
      lies in.  The CEA thermodynamic intervals are 
      [200-1,000], [1,000-6,000], [6,000-20,000] K
     */
    unsigned int interval (const libMesh::Real T) const;

    
    //! @returns a pointer to the coefficients in the interval specified.
    /*!   
      The CEA-style equilibrium curve fits are defined in terms of
      _n_coeffs coefficients for each range fit.
    */
    inline
    const libMesh::Real* coefficients(const unsigned int interval) const
    {
      libmesh_assert_less( interval, this->n_intervals() );
      libmesh_assert_less_equal( _n_coeffs*(interval+1), _coefficients.size() );
      
      return &_coefficients[_n_coeffs*interval];
    }

  protected:

    //! The number of coefficients in each interval
    const unsigned int _n_coeffs;

    //! The coefficient data
    /*!
      The coeffcients are packed in linear ordering. That is,
      a0-a9 for the first interval, a0-a9 for the second interval,
      and so on.
     */
    const std::vector<libMesh::Real> _coefficients;
  };

} // namespace GRINS

#endif //GRINS_CEA_CURVE_FIT_H
