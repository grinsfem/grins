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

// This class
#include "grins/curvilinear_solid_mechanics.h"

// GRINS
#include "grins_config.h"
#include "grins/assembly_context.h"

// libMesh
#include "libmesh/getpot.h"
#include "libmesh/fem_system.h"

namespace GRINS
{
  template<unsigned int Dim>
  CurvilinearSolidMechanics<Dim>::CurvilinearSolidMechanics( const GRINS::PhysicsName & physics_name,
                                                             const GRINS::PhysicsName & core_physics_name,
                                                             const GetPot & input )
    : SolidMechanicsAbstract<Dim>(physics_name,core_physics_name,input)
  {
    static_assert( Dim <= 2, "Error: CurvilinearSolidMechanics defined only for Dim <= 2");
  }

  template<unsigned int Dim>
  void CurvilinearSolidMechanics<Dim>::init_context( AssemblyContext & context )
  {
    this->get_fe(context)->get_JxW();
    this->get_fe(context)->get_phi();
    this->get_fe(context)->get_xyz();
    this->get_fe(context)->get_dphidxi();
    this->get_fe(context)->get_dphideta();

    // Need for constructing metric tensors
    this->get_fe(context)->get_dxyzdxi();
    this->get_fe(context)->get_dxyzdeta();
    this->get_fe(context)->get_dxidx();
    this->get_fe(context)->get_dxidy();
    this->get_fe(context)->get_detadx();
    this->get_fe(context)->get_detady();

    // Only need the z-components if we're embedded in 3-dimensions
    if( Dim == 2 && this->_disp_vars.dim() > 2)
      {
        this->get_fe(context)->get_dxidz();
        this->get_fe(context)->get_detadz();
      }
  }

  template class CurvilinearSolidMechanics<1>;
  template class CurvilinearSolidMechanics<2>;

} // end namespace GRINS
