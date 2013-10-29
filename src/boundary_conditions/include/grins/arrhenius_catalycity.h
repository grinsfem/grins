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

#ifndef GRINS_ARRHENIUS_CATALYCITY_H
#define GRINS_ARRHENIUS_CATALYCITY_H

// GRINS
#include "grins/catalycity_base.h"

namespace GRINS
{
  class ArrheniusCatalycity : public CatalycityBase
  {
  public:

    ArrheniusCatalycity( const libMesh::Real gamma0, const libMesh::Real Ta );

    virtual ~ArrheniusCatalycity();

    virtual libMesh::Real operator()( const libMesh::Real T ) const;

    virtual libMesh::Real dT( const libMesh::Real T ) const;
    
    virtual void set_params( const std::vector<libMesh::Real>& params );

    //! Creates a new copy of the current class.
    /*! A raw pointer is returned and it is assumed the user will take ownership
        and worry about memory management. */
    virtual CatalycityBase* clone() const;
    
  protected:

    libMesh::Real _gamma0;

    libMesh::Real _Ta;

  private:

    ArrheniusCatalycity();

  };

} // end namespace GRINS

#endif // GRINS_ARRHENIUS_CATALYCITY_H
