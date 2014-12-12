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


#ifndef GRINS_SPALART_ALLMARAS_VISCOSITY_H
#define GRINS_SPALART_ALLMARAS_VISCOSITY_H

//GRINS
#include "grins/assembly_context.h"
#include "grins/parameter_user.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/fem_system.h"
#include "libmesh/quadrature.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/function_base.h"

#include "grins/constant_viscosity.h"
#include "grins/parsed_viscosity.h"

class GetPot;

namespace GRINS
{
  template<class Viscosity>
    class SpalartAllmarasViscosity : public ParameterUser
  {
  public:

    SpalartAllmarasViscosity( const GetPot& input );
    ~SpalartAllmarasViscosity();
    
    libMesh::Real operator()(AssemblyContext& context, unsigned int qp) const;

  protected:
    
    //! Viscosity object (so we have access to the physical viscosity)
    Viscosity _mu;
    
  private:

    SpalartAllmarasViscosity();
       
  };

  /* ------------------------- Inline Functions -------------------------*/  
  inline
    template<class Mu>
    libMesh::Real SpalartAllmarasViscosity<Mu>::operator()(AssemblyContext& context, unsigned int qp) const
  {
    // FIXME: We should be getting the variable index to get the qps from the context
    // not hardcode it to be 0
    const std::vector<libMesh::Point>& x = context.get_element_fe(0)->get_xyz();

    const libMesh::Point& x_qp = x[qp];

    // Compute the value of the total viscosity and return it
    libMesh::Number _mu_value = context.interior_value(_turb_vars.turb_var(),x_qp) + this->_mu(context, qp); // Turbulent viscosity + physical viscosity
    
    return _mu_value;
  }    

} // end namespace GRINS

#endif // GRINS_SPALART_ALLMARAS_VISCOSITY_H
