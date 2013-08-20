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

#ifndef GRINS_JOULE_HEATING_H
#define GRINS_JOULE_HEATING_H

// GRINS
#include "grins/physics.h"

//libMesh
#include "libmesh/enum_order.h"
#include "libmesh/enum_fe_family.h"

namespace GRINS
{  
  //! Adds Joule heating source term
  /*!
    This class implements the Joule heating energy term 
   */
  class JouleHeating : public Physics
  {
  public:
    
    JouleHeating( const std::string& physics_name, const GetPot& input );

    ~JouleHeating();

    //! Initialization of AxisymmetricJouleHeating variables
    virtual void init_variables( libMesh::FEMSystem* system );

    //! Source term contribution for AxisymmetricJouleHeating
    /*! This is the main part of the class. This will add the source term to
        the AxisymmetricIncompNavierStokes class.
     */
    virtual void element_time_derivative( bool compute_jacobian,
					  libMesh::FEMContext& context,
					  CachedValues& cache );

  protected:

    //! Element type, read from input
    libMeshEnums::FEFamily _T_FE_family, _V_FE_family;

    //! Temperature element order, read from input
    libMeshEnums::Order _T_order, _V_order;

    //! Name of Temperature
    std::string _T_var_name;

    //! Name of electric potential
    std::string _V_var_name;

    //!
    libMesh::Number _sigma;

    //! Physical dimension of problem
    unsigned int _dim;

    //! Index for Temperature field
    VariableIndex _T_var;

    //! Index for electric potential
    VariableIndex _V_var;

  private:

    JouleHeating();

  };

} // end namespace GRINS

#endif //GRINS_AXISYM_JOULE_HEATING_H
