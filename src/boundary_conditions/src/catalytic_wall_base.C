//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2016 Paul T. Bauman, Roy H. Stogner
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
#include "grins/catalytic_wall_base.h"

// GRINS
#include "grins/math_constants.h"
#include "grins/assembly_context.h"
#include "grins/cached_values.h"

// libMesh
#include "libmesh/fem_system.h"

namespace GRINS
{
  template<typename Chemistry>
  CatalyticWallBase<Chemistry>::CatalyticWallBase( SharedPtr<Chemistry>& chem,
                                                   SharedPtr<CatalycityBase>& gamma,
                                                   const std::vector<VariableIndex>& species_vars,
                                                   VariableIndex T_var,
                                                   libMesh::Real p0,
                                                   unsigned int reactant_species_idx)
    : _chem_ptr(chem),
      _chemistry(*(chem.get())),// This will be removed after NeumannBC refactoring
      _gamma_ptr(gamma),
      _C( std::sqrt( chem->R(reactant_species_idx)/(GRINS::Constants::two_pi) ) ),
      _species_vars(species_vars),
      _T_var(T_var),
      _p0(p0)
  {}

  template<typename Chemistry>
  CatalyticWallBase<Chemistry>::CatalyticWallBase( const Chemistry& chemistry,
                                                   CatalycityBase& gamma,
                                                   const unsigned int reactant_species_idx )
    : _chemistry(chemistry),
      _gamma_s( gamma.clone() ),
      _C( std::sqrt( chemistry.R(reactant_species_idx)/(GRINS::Constants::two_pi) ) )
  {}

  template<typename Chemistry>
  void CatalyticWallBase<Chemistry>::set_catalycity_params( const std::vector<libMesh::Real>& params )
  {
    if(_gamma_s)
      {
        libmesh_deprecated();
        _gamma_s->set_params( params );
      }
    else if( _gamma_ptr)
      _gamma_ptr->set_params( params );
    else
      libmesh_error();
  }

  template<typename Chemistry>
  void CatalyticWallBase<Chemistry>::register_parameter(const std::string & param_name,
                                                        libMesh::ParameterMultiAccessor< libMesh::Number > & param_pointer) const
  {
    _gamma_ptr->register_parameter(param_name,param_pointer);
  }
} // end namespace GRINS
