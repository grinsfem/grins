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
#include "grins/catalytic_wall_neumann_bc_factory_base.h"

// GRINS
#include "grins/catalycity_factory_abstract.h"
#include "grins/gas_recombination_catalytic_wall_neumann_bc_factory_impl.h"
#include "grins/gas_solid_catalytic_wall_neumann_bc_factory_impl.h"

namespace GRINS
{
  template<typename ImplType>
  std::shared_ptr<NeumannBCAbstract>
  CatalyticWallNeumannBCFactoryBase<ImplType>::build_neumann_func( const GetPot& input,
                                                                   MultiphysicsSystem& /*system*/,
                                                                   const FEVariablesBase& fe_var,
                                                                   const std::string& section )
  {
    std::string reaction = this->parse_reaction(input,section);

    // Parse and construct the corresponding catalycity
    std::shared_ptr<CatalycityBase> gamma_ptr = this->build_catalycity( input, section );

    std::string material;
    this->extract_material( fe_var, material );

    // Parse thermodynamic pressure
    /*! \todo We're assuming constant thermodynamic pressure */
    libMesh::Real p0 = this->parse_thermo_pressure(input,material);

    std::string thermochem_lib = this->parse_thermochem_model( input, material );

    return this->build_catalytic_wall_common(input,fe_var,material,reaction,gamma_ptr,p0,thermochem_lib);
  }

  template<typename ImplType>
  std::string CatalyticWallNeumannBCFactoryBase<ImplType>::parse_reaction( const GetPot& input,
                                                                           const std::string& section ) const
  {
    std::string reaction_input_str = section+"/catalytic_reaction";

    // First make sure the input reaction is there
    if(!input.have_variable(reaction_input_str))
      libmesh_error_msg("ERROR: Could not find input for "+reaction_input_str+" !\n");

    // Make sure there's only one
    if( input.vector_variable_size(reaction_input_str) != 1 )
      libmesh_error_msg("ERROR: Can only specify one catalytic_reaction!\n");

    return input( reaction_input_str, std::string("DIE!") );
  }

  template<typename ImplType>
  libMesh::Real CatalyticWallNeumannBCFactoryBase<ImplType>::parse_thermo_pressure( const GetPot& input,
                                                                                    const std::string& material ) const
  {
    std::string thermo_press_input_str =
      "Materials/"+material+"/ThermodynamicPressure/value";

    if( !input.have_variable(thermo_press_input_str) )
      libmesh_error_msg("ERROR: Could not find variable "+thermo_press_input_str+"!");

    return input(thermo_press_input_str, std::numeric_limits<libMesh::Real>::max() );
  }

  template<typename ImplType>
  std::string CatalyticWallNeumannBCFactoryBase<ImplType>::parse_thermochem_model( const GetPot& input,
                                                                                   const std::string& material ) const
  {
    std::string thermochem_input_str = "Materials/"+material+"/GasMixture/thermochemistry_library";

    if( !input.have_variable(thermochem_input_str) )
      libmesh_error_msg("ERROR: Could not find input option "+thermochem_input_str+" !");

    return input(thermochem_input_str, std::string("DIE!") );
  }

  template<typename ImplType>
  std::shared_ptr<CatalycityBase>
  CatalyticWallNeumannBCFactoryBase<ImplType>::build_catalycity( const GetPot& input,
                                                                 const std::string& section ) const
  {
    CatalycityFactoryAbstract::set_getpot(input);
    CatalycityFactoryAbstract::set_section(section);

    std::string catalycity_input_str = section+"/catalycity_type";
    if( !input.have_variable(catalycity_input_str) )
      libmesh_error_msg("ERROR: Could not find input variable "+catalycity_input_str+" !\n");

    std::string catalycity_type = input(catalycity_input_str, "none");

    std::unique_ptr<CatalycityBase> catalycity_ptr = CatalycityFactoryAbstract::build(catalycity_type);

    // We need to return a std::shared_ptr
    return std::shared_ptr<CatalycityBase>( catalycity_ptr.release() );
  }

  template class CatalyticWallNeumannBCFactoryBase<GasRecombinationCatalyticWallNeumannBCFactoryImpl>;
  template class CatalyticWallNeumannBCFactoryBase<GasSolidCatalyticWallNeumannBCFactoryImpl>;

} // end namespace GRINS
