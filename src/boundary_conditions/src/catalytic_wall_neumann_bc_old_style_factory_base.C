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
#include "grins/catalytic_wall_neumann_bc_old_style_factory_base.h"

// GRINS
#include "grins/materials_parsing.h"
#include "grins/catalycity_factory_old_style_base.h"
#include "grins/gas_recombination_catalytic_wall_neumann_bc_factory_impl.h"
#include "grins/gas_solid_catalytic_wall_neumann_bc_factory_impl.h"
#include "grins/string_utils.h"

namespace GRINS
{
  template<typename ImplType>
  std::shared_ptr<NeumannBCAbstract>
  CatalyticWallNeumannBCOldStyleFactoryBase<ImplType>::build_neumann_func( const GetPot& input,
                                                                           MultiphysicsSystem& /*system*/,
                                                                           const FEVariablesBase& fe_var,
                                                                           const std::string& section )
  {
    std::string reaction = this->parse_reaction(input,section);

    // Parse and construct the corresponding catalycity
    std::shared_ptr<CatalycityBase> gamma_ptr = this->build_catalycity( input, section,
                                                                  this->reactant_for_catalycity(reaction) );

    std::string material;
    this->extract_material( fe_var, material );

    // Parse thermodynamic pressure
    /*! \todo We're assuming constant thermodynamic pressure */
    libMesh::Real p0 = this->parse_thermo_pressure(input,material);

    std::string thermochem_lib;
    MaterialsParsing::thermochemistry_lib( input,
                                             PhysicsNaming::reacting_low_mach_navier_stokes(),
                                             thermochem_lib );

    return this->build_catalytic_wall_common(input,fe_var,material,reaction,gamma_ptr,p0,thermochem_lib);
  }

  template<typename ImplType>
  std::string CatalyticWallNeumannBCOldStyleFactoryBase<ImplType>::parse_reaction( const GetPot& input,
                                                                                   const std::string& section ) const
  {
    libmesh_assert_equal_to(NeumannBCFactoryAbstract::_bc_ids->size(), 1 );

    std::string bc_id_string = StringUtilities::T_to_string<BoundaryID>( *(_bc_ids->begin()) );

    std::string prefix_str = this->catalytic_wall_prefix_str();

    std::string reaction_input_str = section+"/"+prefix_str+"_"+bc_id_string;

    // First make sure the input reaction is there
    if(!input.have_variable(reaction_input_str))
      libmesh_error_msg("ERROR: Could not find input for "+reaction_input_str+" !\n");

    // Make sure there's only one
    if( input.vector_variable_size(reaction_input_str) != 1 )
      libmesh_error_msg("ERROR: Can only specify one catalytic_reaction!\n");

    return input( reaction_input_str, std::string("DIE!") );
  }

  template<typename ImplType>
  libMesh::Real CatalyticWallNeumannBCOldStyleFactoryBase<ImplType>::parse_thermo_pressure( const GetPot& input,
                                                                                            const std::string& material ) const
  {
    std::string thermo_press_input_str =
      "Materials/"+material+"/ThermodynamicPressure/value";

    std::string thermo_press_input_str_old_style =
      "Physics/"+PhysicsNaming::reacting_low_mach_navier_stokes()+"/p0";

    if( input.have_variable(thermo_press_input_str) && input.have_variable(thermo_press_input_str_old_style) )
      libmesh_error_msg("ERROR: Cannot specify both "+thermo_press_input_str+" and "+thermo_press_input_str_old_style);

    libMesh::Real invalid_real =  std::numeric_limits<libMesh::Real>::max();
    libMesh::Real p0 = invalid_real;

    if( input.have_variable(thermo_press_input_str) )
      p0 = input(thermo_press_input_str,invalid_real);

    else if( input.have_variable(thermo_press_input_str_old_style) )
      p0 = input(thermo_press_input_str_old_style,invalid_real);

    else
      libmesh_error_msg("ERROR: Could not valid input for thermodynamic pressure!");

    return p0;
  }

  template<typename ImplType>
  std::shared_ptr<CatalycityBase>
  CatalyticWallNeumannBCOldStyleFactoryBase<ImplType>::build_catalycity( const GetPot& input,
                                                                         const std::string& section,
                                                                         const std::string& reactant ) const
  {
    libmesh_assert_equal_to(NeumannBCFactoryAbstract::_bc_ids->size(), 1 );
    std::string bc_id_str = StringUtilities::T_to_string<BoundaryID>( *(_bc_ids->begin()) );

    CatalycityFactoryOldStyleBase::set_getpot(input);
    CatalycityFactoryOldStyleBase::set_section(section);
    CatalycityFactoryOldStyleBase::set_reactant(reactant);
    CatalycityFactoryOldStyleBase::set_bc_id(bc_id_str);

    std::string catalycity_input_str = section+"/gamma_"+reactant+"_"+bc_id_str+"_type";
    if( !input.have_variable(catalycity_input_str) )
      libmesh_error_msg("ERROR: Could not find input variable "+catalycity_input_str+" !\n");

    std::string catalycity_type = input(catalycity_input_str, "none");
    catalycity_type += "_old_style";

    std::unique_ptr<CatalycityBase> catalycity_ptr = CatalycityFactoryOldStyleBase::build(catalycity_type);

    // We need to return a std::shared_ptr
    return std::shared_ptr<CatalycityBase>( catalycity_ptr.release() );
  }

  template class CatalyticWallNeumannBCOldStyleFactoryBase<GasRecombinationCatalyticWallNeumannBCFactoryImpl>;
  template class CatalyticWallNeumannBCOldStyleFactoryBase<GasSolidCatalyticWallNeumannBCFactoryImpl>;

} // end namespace GRINS
