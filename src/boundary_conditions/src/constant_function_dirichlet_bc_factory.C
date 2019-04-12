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
#include "grins/constant_function_dirichlet_bc_factory.h"

// GRINS
#include "grins/string_utils.h"
#include "grins/multicomponent_variable.h"
#include "grins/variable_warehouse.h"
#include "grins/chemistry_builder.h"

#ifdef GRINS_HAVE_CANTERA
#include "grins/cantera_mixture.h"
#endif

#ifdef GRINS_HAVE_ANTIOCH
#include "grins/antioch_chemistry.h"
#endif

// libMesh
#include "libmesh/composite_function.h"
#include "libmesh/const_function.h"
#include "libmesh/zero_function.h"

namespace GRINS
{
  std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
  ConstantFunctionDirichletBCFactory::build_func( const GetPot& input,
                                                  MultiphysicsSystem& system,
                                                  std::vector<std::string>& var_names,
                                                  const std::string& section )
  {
    libmesh_assert( !var_names.empty() );

    //! This is really a "composite" function. We'll cast in the helper functions.
    std::unique_ptr<libMesh::FunctionBase<libMesh::Number> > all_funcs;

    // We're given the active variables in var_names. Let's first check
    // which ones the user actually set in the input file.
    // If there's only one variable in var_names, then check_for_vars will
    // error if the user didn't set it, so we don't need to consider that
    // case here.
    std::set<std::string> vars_found;
    {
      std::vector<std::string> vars_to_search_for(var_names.size());
      this->set_vars_to_search_for(section,var_names,vars_to_search_for);
      this->check_for_vars(input,section,vars_to_search_for,&vars_found);
    }

    all_funcs.reset( this->build_composite_func().release() );

    // Cast to raw CompositeFunction for convenience
    libMesh::CompositeFunction<libMesh::Number>& composite_func =
      libMesh::cast_ref<libMesh::CompositeFunction<libMesh::Number>&>(*(all_funcs.get()));

    std::set<std::string> vars_added;
    this->add_found_vars(input, system, section, vars_found, composite_func, vars_added);

    // Now add all the other vars that weren't added as ZeroFunctions
    for( std::vector<std::string>::const_iterator var = var_names.begin();
         var < var_names.end(); ++var )
      {
        if( vars_added.find(*var) == vars_added.end() )
          {
            std::vector<VariableIndex> var_idx(1,system.variable_number(*var));
            composite_func.attach_subfunction( libMesh::ZeroFunction<libMesh::Number>(),var_idx);
          }
      }

    return all_funcs;
  }

  void ConstantFunctionDirichletBCFactory::add_found_vars( const GetPot& input,
                                                           MultiphysicsSystem& system,
                                                           const std::string& section,
                                                           const std::set<std::string>& vars_found,
                                                           libMesh::CompositeFunction<libMesh::Number>& composite_func,
                                                           std::set<std::string>& vars_added ) const
  {
    libMesh::Number invalid_num = std::numeric_limits<libMesh::Number>::max();

    for( std::set<std::string>::const_iterator var = vars_found.begin();
         var != vars_found.end(); ++var )
      {
        std::vector<VariableIndex> var_idx(1,system.variable_number(*var));

        libMesh::Number value = input(section+"/"+(*var),invalid_num);

        libMesh::ConstFunction<libMesh::Number> const_func(value);
        composite_func.attach_subfunction(const_func,var_idx);
      }

    vars_added = vars_found;
  }

  void MoleFractionsDirichletBCFactory::extract_species_name( const std::string& var_name,
                                                              const std::string& prefix,
                                                              std::string& species_name ) const
  {
    std::vector<std::string> split_name;
    StringUtilities::split_string(var_name,prefix,split_name);

    // split_string won't add the prefix, since it was used as the delimiter, so
    // split_name should just have the lingering species name
    libmesh_assert_equal_to(split_name.size(), 1);

    species_name = split_name[0];
  }

  void MoleFractionsDirichletBCFactory::set_vars_to_search_for( const std::string& section,
                                                                const std::vector<std::string>& var_names,
                                                                std::vector<std::string>&vars_to_search_for ) const
  {
    libmesh_assert_equal_to(var_names.size(),vars_to_search_for.size());

    // Strip out the Variable name from the section
    std::string var_section = extract_var_section(section);

    // This only makes sense for SpeciesMassFractionsVariable in the VariableWarehouse.
    // This call will error out if it's not there.
    const SpeciesMassFractionsVariable& species_fe_var =
      GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>
      (var_section);

    const std::string& prefix = species_fe_var.prefix();
    for( unsigned int v = 0; v < var_names.size(); v++ )
      {
        std::string species_name;
        this->extract_species_name(var_names[v],prefix,species_name);
        vars_to_search_for[v] = "X_"+species_name;
      }
  }

  // To avoid compiler warnings without GRINS or Cantera
#if defined(GRINS_HAVE_ANTIOCH) || defined(GRINS_HAVE_CANTERA)
  void MoleFractionsDirichletBCFactory::add_found_vars( const GetPot& input,
                                                        MultiphysicsSystem& /*system*/,
                                                        const std::string& section,
                                                        const std::set<std::string>& vars_found,
                                                        libMesh::CompositeFunction<libMesh::Number>& composite_func,
                                                        std::set<std::string>& vars_added ) const
#else
  void MoleFractionsDirichletBCFactory::add_found_vars( const GetPot& input,
                                                        MultiphysicsSystem& /*system*/,
                                                        const std::string& section,
                                                        const std::set<std::string>& /*vars_found*/,
                                                        libMesh::CompositeFunction<libMesh::Number>& /*composite_func*/,
                                                        std::set<std::string>& /*vars_added*/ ) const
#endif
  {
    // Strip out the Variable name from the section
    std::string var_section = extract_var_section(section);

    // This only makes sense for SpeciesMassFractionsVariable in the VariableWarehouse.
    // This call will error out if it's not there.
    const SpeciesMassFractionsVariable& species_fe_var =
      GRINSPrivate::VariableWarehouse::get_variable_subclass<SpeciesMassFractionsVariable>
      (var_section);

    const std::string& material = species_fe_var.material();

    std::string thermochem_input_str = "Materials/"+material+"/GasMixture/thermochemistry_library";

    if( !input.have_variable(thermochem_input_str) )
      libmesh_error_msg("ERROR: Could not find input option "+thermochem_input_str+" !");

    const std::string thermochem_lib = input(thermochem_input_str, std::string("DIE!") );

    if( thermochem_lib == "cantera" )
      {
#ifdef GRINS_HAVE_CANTERA
        this->add_mole_frac_to_mass_frac<CanteraMixture>(input,section,vars_found,material,
                                                         species_fe_var,composite_func,vars_added);
#else
        libmesh_error_msg("Error: Cantera not enabled in this configuration. Reconfigure using --with-cantera option.");
#endif
      }
    else if( thermochem_lib == "antioch" )
      {
#ifdef GRINS_HAVE_ANTIOCH
        this->add_mole_frac_to_mass_frac<AntiochChemistry>(input,section,vars_found,material,
                                                           species_fe_var,composite_func,vars_added);
#else
        libmesh_error_msg("Error: Antioch not enabled in this configuration. Reconfigure using --with-antioch option.");
#endif
      }
    else
      libmesh_error_msg("ERROR: Invalid thermochemistry library "+thermochem_lib+"!");
  }

  template<typename ChemistryType>
  void MoleFractionsDirichletBCFactory::add_mole_frac_to_mass_frac(const GetPot& input,
                                                                   const std::string& section,
                                                                   const std::set<std::string>& vars_found,
                                                                   const std::string& material,
                                                                   const SpeciesMassFractionsVariable& species_fe_var,
                                                                   libMesh::CompositeFunction<libMesh::Number>& composite_func,
                                                                   std::set<std::string>& vars_added) const
  {
    unsigned int n_vars_found = vars_found.size();

    // Parse in all the species mole fracs that are in the input (it is assumed non-specified are 0)
    std::vector<libMesh::Number> species_mole_fracs(n_vars_found);
    libMesh::Number invalid_num = std::numeric_limits<libMesh::Number>::max();
    {
      unsigned int count = 0;
      for(std::set<std::string>::const_iterator var = vars_found.begin();
          var != vars_found.end(); ++var )
        {
          species_mole_fracs[count] = input(section+"/"+(*var),invalid_num);
          count++;
        }
    }

    // Make sure mole fracs sum to 1
    libMesh::Number sum = 0.0;
    for(unsigned int v = 0; v < n_vars_found; v++ )
      sum += species_mole_fracs[v];

    libMesh::Number tol = std::numeric_limits<libMesh::Number>::epsilon()*10;
    if( std::abs(sum-1.0) > tol )
      libmesh_error_msg("ERROR: Mole fractions do not sum to 1! Found sum = "+StringUtilities::T_to_string<libMesh::Number>(sum));

    // Extract species names
    std::vector<std::string> species_names(n_vars_found);
    {
      unsigned int count = 0;
      for(std::set<std::string>::const_iterator var = vars_found.begin();
          var != vars_found.end(); ++var )
        {
          std::vector<std::string> split_name;
          // vars_found should have the form "X_<species name>"
          StringUtilities::split_string((*var),"_",split_name);
          libmesh_assert_equal_to(split_name[0],std::string("X"));
          libmesh_assert_equal_to(split_name.size(),2);
          species_names[count] = split_name[1];
          count++;
        }
    }

    // Now convert to mass frac and add to composite function
    /*! \todo We should have a ChemsitryWarehouse or something to just grab this from one place
      instead of rebuilding. */
    ChemistryBuilder chem_builder;
    std::unique_ptr<ChemistryType> chem_ptr;
    chem_builder.build_chemistry(input,material,chem_ptr);

    const ChemistryType & chem = *chem_ptr;

    libMesh::Real M = 0.0;
    for(unsigned int v = 0; v < n_vars_found; v++ )
      {
        unsigned int s = chem.species_index(species_names[v]);
        M += species_mole_fracs[v]*chem.M(s);
      }

    const std::string& prefix = species_fe_var.prefix();

    for(unsigned int v = 0; v < n_vars_found; v++ )
      {
        // Finish computing species mass fraction
        unsigned int s = chem.species_index(species_names[v]);
        libMesh::Number species_mass_fracs = species_mole_fracs[v]*chem.M(s)/M;

        // Add the function
        std::vector<VariableIndex> var_idx(1,species_fe_var.species(s));
        libMesh::ConstFunction<libMesh::Number> const_func(species_mass_fracs);
        composite_func.attach_subfunction(const_func,var_idx);

        // Log that we added this variable
        vars_added.insert(prefix+species_names[v]);
      }
  }

  std::string MoleFractionsDirichletBCFactory::extract_var_section( const std::string& section ) const
  {
    std::vector<std::string> tokens;
    StringUtilities::split_string(section,"/",tokens);
    return tokens.back();
  }

  // Instantiate
#ifdef GRINS_HAVE_CANTERA
  template void MoleFractionsDirichletBCFactory::add_mole_frac_to_mass_frac<CanteraMixture>(const GetPot&,const std::string&,const std::set<std::string>&,const std::string&,const SpeciesMassFractionsVariable&,libMesh::CompositeFunction<libMesh::Number>&,std::set<std::string>& ) const;
#endif

#ifdef GRINS_HAVE_ANTIOCH
  template void MoleFractionsDirichletBCFactory::add_mole_frac_to_mass_frac<AntiochChemistry>(const GetPot&,const std::string&,const std::set<std::string>&,const std::string&,const SpeciesMassFractionsVariable&,libMesh::CompositeFunction<libMesh::Number>&,std::set<std::string>& ) const;
#endif

} // end namespace GRINS
