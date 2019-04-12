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

#ifndef GRINS_CONSTANT_FUNCTION_DIRICHLET_BC_FACTORY_H
#define GRINS_CONSTANT_FUNCTION_DIRICHLET_BC_FACTORY_H

// GRINS
#include "grins/dirichlet_bc_factory_function_base.h"
#include "grins/parsed_function_factory_helper.h"

// libMesh
#include "libmesh/function_base.h"

namespace GRINS
{
  // Foward declarations
  class SpeciesMassFractionsVariable;

  //! Constructs ConstFunction objects for Dirichlet boundary conditions
  class ConstantFunctionDirichletBCFactory : public DirichletBCFactoryFunctionBase<libMesh::FunctionBase<libMesh::Number> >,
                                             public ParsedFunctionFactoryHelper<libMesh::FunctionBase<libMesh::Number> >
  {
  public:

    ConstantFunctionDirichletBCFactory( const std::string& bc_type_name )
      : DirichletBCFactoryFunctionBase<libMesh::FunctionBase<libMesh::Number> >(bc_type_name),
      ParsedFunctionFactoryHelper<libMesh::FunctionBase<libMesh::Number> >()
    {}

    ~ConstantFunctionDirichletBCFactory(){};

  protected:
    //! Builds ConstantFunction objects for boundary conditions
    /*! The variable names passed in will correspond to only a single
      VariableBase object, e.g. Velocity. The expected behavior is
      that if the user didn't specify a value for all the variables,
      then the unspecified variables will be set to zero. However,
      the user must've set at least one. The section arguments
      corresponds to the section to parse for the variables in the
      input file, e.g. input(section+"/"+var_names[0]). */
    virtual std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
    build_func( const GetPot& input,
                MultiphysicsSystem& system,
                std::vector<std::string>& var_names,
                const std::string& section );

    //! Adds the vars_found to the composite_func
    /*! It is expected that this function returns the names of the vars actually added
      in vars_added. While the default implementation will have vars_found = vars_added,
      subclasses may actually search for different variables and use them to add the actual
      system variables. */
    virtual void add_found_vars(const GetPot& input,
                                MultiphysicsSystem& system,
                                const std::string& section,
                                const std::set<std::string>& vars_found,
                                libMesh::CompositeFunction<libMesh::Number>& composite_func,
                                std::set<std::string>& vars_added) const;

    //! Set the vars_to_search_for, based on var_names
    /*! In the default implementation, vars_to_search_for = var_names, but subclasses
      may override this behavior. If this function is overridden, then add_found_vars
      will also need to be overridden. You may assume vars_to_search_for has been
      sized to match var_names. */
    virtual void set_vars_to_search_for( const std::string& /*section*/,
                                         const std::vector<std::string>& var_names,
                                         std::vector<std::string>& vars_to_search_for ) const
    { libmesh_assert_equal_to(var_names.size(),vars_to_search_for.size());
      vars_to_search_for = var_names; }

  };

  //! Parses mole fraction values and converts to mass fractions
  /*! Only valid for SpeciesMassFraction FEVariables. */
  class MoleFractionsDirichletBCFactory : public ConstantFunctionDirichletBCFactory
  {
  public:
    MoleFractionsDirichletBCFactory( const std::string& bc_type_name )
      : ConstantFunctionDirichletBCFactory(bc_type_name)
    {}

    ~MoleFractionsDirichletBCFactory(){}

  protected:

    //! Here, we're expected vars_found to correspond to mole fractions and we'll add mass fractions
    /*! vars_found should have things like X_N, etc. The prefix will be "X_" for mole fractions.
      Then, we'll add mass fractions, using the corresponding names in the SpeciesMassFractionsVariable.
      We'll match them based on the species names. */
    virtual void add_found_vars(const GetPot& input,
                                MultiphysicsSystem& system,
                                const std::string& section,
                                const std::set<std::string>& vars_found,
                                libMesh::CompositeFunction<libMesh::Number>& composite_func,
                                std::set<std::string>& vars_added) const;

    //! We'll search for mole fractions: X_<species_name>.
    /*! We extract the species name from var_names and then set the corresponding
      vars_to_search_for to X_<species_name>. */
    virtual void set_vars_to_search_for( const std::string& section,
                                         const std::vector<std::string>& var_names,
                                         std::vector<std::string>&vars_to_search_for ) const;

    template<typename ChemistryType>
    void add_mole_frac_to_mass_frac(const GetPot& input,
                                    const std::string& section,
                                    const std::set<std::string>& vars_found,
                                    const std::string& material,
                                    const SpeciesMassFractionsVariable& species_fe_var,
                                    libMesh::CompositeFunction<libMesh::Number>& composite_func,
                                    std::set<std::string>& vars_added) const;

    void extract_species_name( const std::string& var_name,
                               const std::string& prefix,
                               std::string& species_name ) const;

    std::string extract_var_section( const std::string& section ) const;

  };

} // end namespace GRINS

#endif // GRINS_CONSTANT_FUNCTION_DIRICHLET_BC_FACTORY_H
