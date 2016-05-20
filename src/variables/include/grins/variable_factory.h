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

#ifndef GRINS_VARIABLE_FACTORY_H
#define GRINS_VARIABLE_FACTORY_H

// GRINS
#include "grins/factory_with_getpot.h"
#include "grins/fe_variables_base.h"

// libMesh
#include "libmesh/factory.h"

namespace GRINS
{
  class VariableFactoryAbstract : public FactoryWithGetPot<FEVariablesBase>
  {
  public:

    VariableFactoryAbstract( const std::string& name )
      : FactoryWithGetPot<FEVariablesBase>(name)
    {}

    ~VariableFactoryAbstract(){};

    virtual libMesh::UniquePtr<FEVariablesBase> create();

    //! Build the variable names for the FEVariablesBase type (name), returned in the std::vector
    /*! Similarly to build(), this will grab the factory subclass and then
        internally call parse_var_names. set_getpot() must called before
        calling this method. */
    static std::vector<std::string> build_var_names( const std::string& name );

    //! Set the variable names before calling create()
    static void set_var_names( const std::vector<std::string>& var_names )
    { _var_names = &var_names; }

    //! Set the variable indices before calling create()
    static void set_var_indices( const std::vector<VariableIndex>& var_indices )
    { _var_indices = &var_indices; }

    //! Set the section for the input file before calling build_var_names()
    /*! It is expected that user just need to call input(var_section+"/<user_option>") */
    static void set_var_section( const std::string& var_section )
    { _var_section = var_section; }

  protected:

    //! Helper function to check required data is set when calling create()
    virtual void check_create_state() const;

    //! Helper function to reset data before next call to create()
    virtual void reset_create_state();

    //! Subclasses implement construction of the FEVariablesBase object using the var_names and var_indices
    /*! This function will be called from within create(), which called from
        VariableFactoryAbstract::build. Note the var_names can be built a priori
        using the VariableFactoryAbstract::build_var_names() method. */
    virtual libMesh::UniquePtr<FEVariablesBase> build_fe_var( const std::vector<std::string>& var_names,
                                                              const std::vector<VariableIndex>& var_indices ) =0;

    // Subclasses implement this to parse the variable component name(s)
    virtual std::vector<std::string> parse_var_names( const GetPot& input, const std::string& var_section ) =0;

    //! Variable component names needed for FEVariableBase construction
    static const std::vector<std::string>* _var_names;

    //! Variable component indices needed for FEVariableBase construction
    static const std::vector<VariableIndex>* _var_indices;

    //! Section of input to parse variable names in build_var_names
    static std::string _var_section;

  private:

    VariableFactoryAbstract();

  };


  //! Factory to build "standard" FEVariablesBase classes
  template<typename VariableType>
  class VariableFactoryBasic : public VariableFactoryAbstract
  {
  public:

    VariableFactoryBasic( const std::string& name )
      : VariableFactoryAbstract(name)
    {}

    ~VariableFactoryBasic(){}

  protected:

    virtual libMesh::UniquePtr<FEVariablesBase> build_fe_var( const std::vector<std::string>& var_names,
                                                              const std::vector<VariableIndex>& var_indices )
    { return libMesh::UniquePtr<FEVariablesBase>( new VariableType(var_names,var_indices) ); }

    //! The basic factory implementation looks in [Variables/<VariableName>/names].
    virtual std::vector<std::string> parse_var_names( const GetPot& input, const std::string& var_section );

  };


  //! Factory to build FEVariablesBase classes that use species names as variables
  /*! Thus, we need a special way to parse the input to figure out what all
      the species names are. */
  template<typename VariableType>
  class SpeciesVariableFactory : public VariableFactoryAbstract
  {
  public:

    SpeciesVariableFactory( const std::string& name )
      : VariableFactoryAbstract(name)
    {}

    ~SpeciesVariableFactory(){}

  protected:

    //! Implementation species variable name parsing
    /*! First, we look for the [Variables/<VariableType>/prefix], which is the prefix for all
        the species variable component names. The, we need to look up the material to figure out
        where to grab the species from. With the material name, then we look up the species names
        and accordingly build up the variable names. */
    virtual std::vector<std::string> parse_var_names( const GetPot& input, const std::string& var_section );

    virtual libMesh::UniquePtr<FEVariablesBase> build_fe_var( const std::vector<std::string>& var_names,
                                                              const std::vector<VariableIndex>& var_indices )
    { return libMesh::UniquePtr<FEVariablesBase>( new VariableType(var_names,var_indices,_prefix,_material) ); }

    std::string _prefix;

    std::string _material;
  };

} // end namespace GRINS

#endif // GRINS_VARIABLE_FACTORY_H
