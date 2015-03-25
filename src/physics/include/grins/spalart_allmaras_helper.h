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


#ifndef GRINS_SPALART_ALLMARAS_HELPER_H
#define GRINS_SPALART_ALLMARAS_HELPER_H

//GRINS
#include "grins/physics.h"
#include "primitive_flow_variables.h"

//Utils
#include "grins/distance_function.h"

namespace GRINS
{
  class SpalartAllmarasHelper
  {
  public:

    SpalartAllmarasHelper(const GetPot& input);

    ~SpalartAllmarasHelper();

    void init_variables( libMesh::FEMSystem* system );
            
    // The vorticity function
    libMesh::Real _vorticity(AssemblyContext& context, unsigned int qp) const;
    
    // The source function \tilde{S}
    libMesh::Real _source_fn( libMesh::Number nu, libMesh::Real mu, libMesh::Real wall_distance, libMesh::Real _vorticity_value) const;

    // The destruction function f_w(nu)
    libMesh::Real _destruction_fn(libMesh::Number nu, libMesh::Real wall_distance, libMesh::Real _S_tilde) const;
    
    //! Spalart Allmaras model constants
    libMesh::Number _cb1, _sigma, _cb2, _cw1;

    //! Constants specific to the calculation of the source function
    libMesh::Number _kappa, _cv1, _cv2, _cv3;

    //! Constants specific to the calculation of the destruction function
    libMesh::Number _r_lin, _c_w2, _c_w3;    

    //! Constants specific to the calculation of the trip function (but used in
    // the source and destruction term)
    libMesh::Number _c_t3, _c_t4;

    //! Constants specific to the calculation of the negative S-A model
    libMesh::Number _c_n1;
    
  protected:

    // Physical dimension
    unsigned int _dim;

    // The flow variables
    PrimitiveFlowVariables _flow_vars;
    
  private:
    SpalartAllmarasHelper();

  };

} //End namespace block

#endif // GRINS_SPALART_ALLMARAS_HELPER_H
