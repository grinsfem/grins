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
#ifndef QOI_FACTORY_H
#define QOI_FACTORY_H

//libMesh
#include "libmesh/getpot.h"

//GRINS
#include "grins/grins_physics_names.h"
#include "grins/qoi_names.h"
#include "grins/qoi_base.h"
#include "grins/average_nusselt_number.h"
#include "grins/vorticity.h"

// shared_ptr
#include "boost/tr1/memory.hpp"

namespace GRINS
{
  class QoIFactory
  {
  public:

    QoIFactory();
    
    virtual ~QoIFactory();

    virtual std::tr1::shared_ptr<QoIBase> build(const GetPot& input);

  protected:

    virtual void add_qoi( const GetPot& input, const std::string& qoi_name, std::tr1::shared_ptr<QoIBase>& qoi );

    virtual void check_qoi_physics_consistency( const GetPot& input,
						const std::string& qoi_name );

    virtual void echo_qoi_list( const std::string& qoi_name );

    void consistency_helper( const std::set<std::string>& requested_physics,
			     const std::set<std::string>& required_physics, 
			     const std::string& qoi_name );

    void consistency_error_msg( const std::string& qoi_name, const std::set<std::string>& required_physics );

  };
}
#endif // QOI_FACTORY_H
