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


#ifndef GRINS_SPALART_ALLMARAS_VISCOSITY_H
#define GRINS_SPALART_ALLMARAS_VISCOSITY_H

//GRINS
#include "grins/parameter_user.h"
#include "grins/spalart_allmaras_parameters.h"

// libMesh
#include "libmesh/libmesh_common.h"
#include "libmesh/point.h"

// libMesh forward declarations
class GetPot;
namespace libMesh
{
  class FEMSystem;
}

namespace GRINS
{
  // Forward declarations
  class AssemblyContext;
  class TurbulenceFEVariables;

  template<class Viscosity>
  class SpalartAllmarasViscosity : public ParameterUser
  {
  public:

    //! Constructor with specified material
    SpalartAllmarasViscosity( const GetPot& input, const std::string& material );

    //! Deprecated constructor
    SpalartAllmarasViscosity( const GetPot& input );
    ~SpalartAllmarasViscosity(){};

    libMesh::Real operator()(AssemblyContext& context, unsigned int qp) const;

    libMesh::Real operator()( const libMesh::Point& p, const libMesh::Real time=0 )
    { return _mu(p,time); }

    // Registers all parameters in this physics and in its property
    // classes
    virtual void register_parameter
    ( const std::string & param_name,
      libMesh::ParameterMultiAccessor<libMesh::Number> & param_pointer )
      const;

  protected:

    //! Viscosity object (so we have access to the physical viscosity)
    Viscosity _mu;

    // These are defined for each physics
    TurbulenceFEVariables& _turbulence_vars;

    SpalartAllmarasParameters _sa_params;

  private:

    SpalartAllmarasViscosity();

  };

} // end namespace GRINS

#endif // GRINS_SPALART_ALLMARAS_VISCOSITY_H
