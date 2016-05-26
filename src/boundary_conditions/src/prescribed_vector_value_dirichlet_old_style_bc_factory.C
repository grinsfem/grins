//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2015 Paul T. Bauman, Roy H. Stogner
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
#include "grins/prescribed_vector_value_dirichlet_old_style_bc_factory.h"

// GRINS
#include "grins/string_utils.h"
#include "grins/species_mass_fracs_fe_variables.h"
#include "grins/variable_warehouse.h"
#include "grins/physics_naming.h"
#include "grins/physics_factory_helper.h"

#ifdef GRINS_HAVE_CANTERA
#include "grins/cantera_mixture.h"
#endif

#ifdef GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

// libMesh
#include "libmesh/composite_function.h"
#include "libmesh/const_function.h"

namespace GRINS
{
  libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Number> >
  PrescribedVectorValueDirichletOldStyleBCFactory::build_func( const GetPot& input,
                                                               MultiphysicsSystem& system,
                                                               std::vector<std::string>& var_names,
                                                               const std::string& section )
  {
    libmesh_assert_equal_to(DirichletBCFactoryAbstract::_bc_ids->size(), 1 );

    std::string bc_id_string = StringUtilities::T_to_string<BoundaryID>( *(_bc_ids->begin()) );

    std::string var_input_string = this->var_input_string();

    std::string input_string = section+"/"+var_input_string+"_"+bc_id_string;

    unsigned int n_comps = input.vector_variable_size(input_string);

    if( var_names.size() > n_comps )
      libmesh_error_msg("ERROR: Insufficient number of variable components in "+input_string+"!");

    libMesh::UniquePtr<libMesh::CompositeFunction<libMesh::Number> >
      remapped_func( new libMesh::CompositeFunction<libMesh::Number> );

    this->add_funcs(input,system,input_string,var_names,*(remapped_func.get()));

    return libMesh::UniquePtr<libMesh::FunctionBase<libMesh::Number> >(remapped_func.release());
  }

  void PrescribedVectorValueDirichletOldStyleBCFactory::add_funcs( const GetPot& input,
                                                                   MultiphysicsSystem& system,
                                                                   const std::string& input_string,
                                                                   const std::vector<std::string>& var_names,
                                                                   libMesh::CompositeFunction<libMesh::Number>& composite_func ) const
  {
    for( unsigned int n = 0; n < var_names.size(); n++ )
      {
        std::vector<VariableIndex> dbc_vars(1,system.variable_number(var_names[n]));
        libMesh::Number value = input(input_string, 0.0, n);
        libMesh::ConstFunction<libMesh::Number> const_func(value);
        composite_func.attach_subfunction(const_func, dbc_vars);
      }
  }

  void PrescribedMoleFractionsDirichletOldStyleBCFactory::add_funcs( const GetPot& input,
                                                                     MultiphysicsSystem& /*system*/,
                                                                     const std::string& input_string,
                                                                     const std::vector<std::string>& var_names,
                                                                     libMesh::CompositeFunction<libMesh::Number>& composite_func ) const
  {
    const unsigned int n_vars = var_names.size();

    // Parse in all the species mole fracs that are in the input
    std::vector<libMesh::Number> species_mole_fracs(n_vars);
    libMesh::Number invalid_num = std::numeric_limits<libMesh::Number>::max();

    if( input.vector_variable_size(input_string) != n_vars )
      libmesh_error_msg("ERROR: Expected "+StringUtilities::T_to_string<unsigned int>(n_vars)+" components in "+input_string);

    for(unsigned int v = 0; v < n_vars; v++ )
      species_mole_fracs[v] = input(input_string,invalid_num,v);

    // Make sure mole fracs sum to 1
    libMesh::Number sum = 0.0;
    for(unsigned int v = 0; v < n_vars; v++ )
      sum += species_mole_fracs[v];

    libMesh::Number tol = std::numeric_limits<libMesh::Number>::epsilon()*10;
    if( std::abs(sum-1.0) > tol )
      libmesh_error_msg("ERROR: Mole fractions do not sum to 1! Found sum = "+StringUtilities::T_to_string<libMesh::Number>(sum));


    // This only makes sense for SpeciesMassFractionsFEVariables in the
    // VariableWarehouse. This call will error out if it's not there.
    const SpeciesMassFractionsFEVariables& species_fe_var =
      GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsFEVariables>
      (VariablesParsing::species_mass_fractions_section());

    std::string thermochem_lib;
    PhysicsFactoryHelper::parse_thermochemistry_model( input,
                                                       PhysicsNaming::reacting_low_mach_navier_stokes(),
                                                       thermochem_lib );

    if( thermochem_lib == "cantera" )
      {
#ifdef GRINS_HAVE_CANTERA
        this->convert_mole_fracs_and_add_to_func<CanteraMixture>(input,
                                                                 species_mole_fracs,
                                                                 species_fe_var,
                                                                 composite_func);
#else
        libmesh_error_msg("Error: Cantera not enabled in this configuration. Reconfigure using --with-cantera option.");
#endif
      }
    else if( thermochem_lib == "antioch" )
      {
#ifdef GRINS_HAVE_ANTIOCH
        this->convert_mole_fracs_and_add_to_func<AntiochChemistry>(input,
                                                                   species_mole_fracs,
                                                                   species_fe_var,
                                                                   composite_func);
#else
        libmesh_error_msg("Error: Antioch not enabled in this configuration. Reconfigure using --with-antioch option.");
#endif
      }
    else
      libmesh_error_msg("ERROR: Invalid thermochemistry library "+thermochem_lib+"!");
  }

  template<typename ChemistryType>
  void PrescribedMoleFractionsDirichletOldStyleBCFactory::
  convert_mole_fracs_and_add_to_func
  (const GetPot& input, const std::vector<libMesh::Number>& species_mole_fracs,
   const SpeciesMassFractionsFEVariables& species_fe_var,
   libMesh::CompositeFunction<libMesh::Number>& composite_func) const
  {
    const std::string& material = species_fe_var.material();

    /*! \todo We should have a ChemsitryWarehouse or something to just
              grab this from one place instead of rebuilding. */
    ChemistryType chem(input,material);

    const unsigned int n_vars = species_mole_fracs.size();
    // Compute M
    libMesh::Real M = 0.0;
    for( unsigned int s = 0; s < n_vars; s++ )
      M += species_mole_fracs[s]*chem.M(s);

    // Convert mole fractions to mass fractions and add to function
    for( unsigned int s = 0; s < n_vars; s++ )
      {
        libMesh::Number species_mass_fracs =species_mole_fracs[s]*chem.M(s)/M;

        std::vector<VariableIndex> var_idx(1,species_fe_var.species(s));

        libMesh::ConstFunction<libMesh::Number> const_func(species_mass_fracs);
        composite_func.attach_subfunction(const_func,var_idx);
      }
  }

  //Instantiate
#ifdef GRINS_HAVE_CANTERA
  template void PrescribedMoleFractionsDirichletOldStyleBCFactory::convert_mole_fracs_and_add_to_func<CanteraMixture>(const GetPot&, const std::vector<libMesh::Number>&, const SpeciesMassFractionsFEVariables&,libMesh::CompositeFunction<libMesh::Number>& ) const;
#endif

#ifdef GRINS_HAVE_ANTIOCH
  template void PrescribedMoleFractionsDirichletOldStyleBCFactory::convert_mole_fracs_and_add_to_func<AntiochChemistry>(const GetPot&, const std::vector<libMesh::Number>&, const SpeciesMassFractionsFEVariables&,libMesh::CompositeFunction<libMesh::Number>& ) const;
#endif

  // Register factories
  PrescribedVelDirichletOldStyleBCFactory grins_factory_prescribed_vel_old_style("prescribed_vel_old_style");
  PrescribedDispDirichletOldStyleBCFactory grins_factory_constant_displacement_old_style("constant_displacement_old_style");
  PrescribedSpeciesDirichletOldStyleBCFactory grins_factory_prescribed_species_old_style("prescribed_species_old_style");
  PrescribedMoleFractionsDirichletOldStyleBCFactory grins_factory_prescribed_mole_fraccs_old_style("prescribed_mole_fracs_old_style");

} // end namespace GRINS
