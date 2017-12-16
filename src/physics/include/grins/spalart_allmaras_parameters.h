//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_SPALART_ALLMARAS_PARAMETERS_H
#define GRINS_SPALART_ALLMARAS_PARAMETERS_H

// GRINS
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/libmesh.h"
class GetPot;

namespace GRINS
{
  //! Encapsulate Spalart-Allmaras model parameters
  /*! This is mostly a container class, but there are a few helper functions
    here that are used in different places in SpalartAllmaras classes. */
  class SpalartAllmarasParameters : public ParameterUser
  {
  public:

    SpalartAllmarasParameters(const GetPot& input);

    ~SpalartAllmarasParameters(){};

    // The source function \f$ \tilde{S} \f$
    libMesh::Real source_fn( libMesh::Number nu, libMesh::Real mu,
                             libMesh::Real wall_distance, libMesh::Real vorticity_value, bool infinite_distance) const;

    // The destruction function \f$ f_w(\nu) \f$
    libMesh::Real destruction_fn( libMesh::Number nu, libMesh::Real wall_distance,
                                  libMesh::Real S_tilde, bool infinite_distance) const;

    //! Helper function
    /*! This expression appears in a couple of places so we provide a function for it*/
    libMesh::Real fv1( libMesh::Real chi ) const;

    libMesh::Real get_kappa() const
    { return _kappa;}

    libMesh::Real get_cv1() const
    { return _cv1;}

    libMesh::Real get_cv2() const
    { return _cv2;}

    libMesh::Real get_cv3() const
    { return _cv3;}

    libMesh::Real get_cb1() const
    { return _cb1;}

    libMesh::Real get_cb2() const
    { return _cb2;}

    libMesh::Real get_sigma() const
    { return _sigma;}

    libMesh::Real get_c_w2() const
    { return _c_w2;}

    libMesh::Real get_c_w3() const
    { return _c_w3;}

    libMesh::Real get_r_lin() const
    { return _r_lin;}

    libMesh::Real get_c_t3() const
    { return _c_t3;}

    libMesh::Real get_c_t4() const
    { return _c_t4;}

    libMesh::Real get_c_n1() const
    { return _c_n1;}

  protected:

    //! Constants specific to the calculation of the source function
    libMesh::Real _kappa, _cv1, _cv2, _cv3;

    //! Spalart Allmaras model constants, the constant _cw1 are calculated, not cached
    libMesh::Real _cb1, _sigma, _cb2;

    //! Constants specific to the calculation of the destruction function
    libMesh::Real _r_lin, _c_w2, _c_w3;

    //! Constants specific to the calculation of the trip function (but used in
    // the source and destruction term)
    libMesh::Real _c_t3, _c_t4;

    //! Constants specific to the calculation of the negative S-A model
    libMesh::Real _c_n1;

  private:

    SpalartAllmarasParameters();

  };

  inline
  libMesh::Real SpalartAllmarasParameters::fv1( libMesh::Real chi ) const
  {
    libMesh::Real chi3 = chi*chi*chi;
    libMesh::Real cv1 = this->get_cv1();
    libMesh::Real cv13 = cv1*cv1*cv1;

    return chi3/(chi3 + cv13);
  }

} // end namespace GRINS

#endif // GRINS_SPALART_ALLMARAS_PARAMETERS_H
