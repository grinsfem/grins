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
#include "grins/neumann_bc_factory_abstract.h"

// GRINS
#include "grins/string_utils.h"

namespace GRINS
{
  std::unique_ptr<NeumannBCContainer> NeumannBCFactoryAbstract::create()
  {
    std::unique_ptr<NeumannBCContainer> container;

    // Only build the container if the Neumann BC is *not* homogeneous
    // Otherwise, we just return and empty (NULL) UniquePtr.
    if(!_is_homogeneous)
      {
        // Make sure all necessary state has been setup
        this->check_state();

        std::shared_ptr<NeumannBCAbstract>
          func = this->build_neumann_func( *_input, *_system, *_fe_var, this->_section );

        libmesh_assert(func);

        container.reset(new NeumannBCContainer(*_bc_ids,*_fe_var,func));

        // Reset state for error checking during next construction
        this->reset_state();
      }

    return container;
  }

  void NeumannBCFactoryAbstract::check_for_flux( const GetPot& input, const std::string& flux_input,
                                                 const std::vector<std::string>& var_names )
  {
    if( !input.have_variable(flux_input) )
      libmesh_error_msg("ERROR: Could not find input specification for "+flux_input+"!");

    unsigned int flux_size = input.vector_variable_size(flux_input);
    if( flux_size != var_names.size() )
      {
        std::string error_msg = "ERROR: Mismatch in size between flux input and variables size!\n";
        error_msg += "       Found flux size      = "+StringUtilities::T_to_string<unsigned int>(flux_size)+"\n";
        error_msg += "       Found variables size = "+StringUtilities::T_to_string<unsigned int>(var_names.size())+"\n";
        libmesh_error_msg(error_msg);
      }
  }

  // Full specialization for the Factory<NeumannBCContainer>
  template<>
  std::map<std::string, FactoryAbstract<NeumannBCContainer>*>&
  FactoryAbstract<NeumannBCContainer>::factory_map()
  {
    static std::map<std::string, FactoryAbstract<NeumannBCContainer>*> _map;
    return _map;
  }

  // Definition of static members
  template<>
  const GetPot* FactoryWithGetPot<NeumannBCContainer>::_input = NULL;

  template<>
  MultiphysicsSystem* BCFactoryAbstract<NeumannBCContainer>::_system = NULL;

  template<>
  const std::set<BoundaryID>* BCFactoryAbstract<NeumannBCContainer>::_bc_ids = NULL;

  template<>
  const FEVariablesBase* BCFactoryAbstract<NeumannBCContainer>::_fe_var = NULL;

  template<>
  std::string BCFactoryAbstract<NeumannBCContainer>::_section = std::string("DIE!");

} // end namespace GRINS
