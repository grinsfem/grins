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

#ifndef GRINS_BC_FACTORY_ABSTRACT_H
#define GRINS_BC_FACTORY_ABSTRACT_H

// GRINS
#include "grins/factory_with_getpot.h"
#include "grins/var_typedefs.h"
#include "grins/fe_variables_base.h"

// libMesh
#include "libmesh/libmesh.h"
#include "libmesh/boundary_info.h" // invalid_id
#include "grins/multiphysics_sys.h"

namespace GRINS
{
  // Forward declarations
  class MultiphysicsSystem;

  template<typename Base>
  class BCFactoryAbstract : public FactoryWithGetPot<Base>
  {
  public:
    BCFactoryAbstract( const std::string& bc_type_name )
      : FactoryWithGetPot<Base>(bc_type_name)
    {}

    ~BCFactoryAbstract(){};

    static void set_system( MultiphysicsSystem& system )
    { _system = &system; }

    //! Boundary id for the current boundary condition section
    static void set_bc_ids( const std::set<BoundaryID>& bc_ids )
    { _bc_ids = &bc_ids; }

    //! Active variable for the current boundary condition
    static void set_fe_var( const FEVariablesBase& fe_var )
    { _fe_var = &fe_var; }

    //! Sets the current section of the input file
    /*! The section here corresponds to the section in the GetPot
      input file, e.g. if we're wanting to parse functions from
      [Variables/SideWall/Velocity/u], then section will be
      "Variables/SideWall/Velocity". */
    static void set_section( const std::string& section )
    { _section = section; }

    static bool have_bc_type( const std::string& bc_type );

  protected:

    //! Ensure that there is only one expression for the [Section/var_name] variable
    /*! When parsing expressions to give to Parsed(FEM)Function, this function checks
      that there's only one expression for the given [Section/var_name] in the input
      file. This also helps protect against white space in the expression, which we
      don't currently support. */
    void check_for_multiple_expressions( const GetPot& input,const std::string& section,
                                         const std::string& var_name ) const;

    void build_var_indices( const MultiphysicsSystem& system,
                            const std::vector<std::string>& var_names,
                            std::vector<VariableIndex>& var_indices ) const;

    //! Helper function to reduce code duplication
    virtual void check_state() const;

    //! Helper function to redue code duplication
    virtual void reset_state();

    /*! We store only a raw pointer here because we *can't* make a copy.
      Otherwise, bad things will happen. We are not taking
      ownership of this, so we need to *not* delete this.*/
    static MultiphysicsSystem* _system;

    //! BoundaryID for constructing a particular boundary condition
    static const std::set<BoundaryID>* _bc_ids;

    //! The FEVariablesBase class associated with the boundary condition being built
    /*! We only build one boundary condition at a time, so this pointer will
      change with each boundary condition construction. We store a raw pointer
      here because we don't own this. Do not delete!*/
    static const FEVariablesBase* _fe_var;

    static std::string _section;

  };

  template<typename Base>
  inline
  void BCFactoryAbstract<Base>::check_for_multiple_expressions( const GetPot& input,const std::string& section,
                                                                const std::string& var_name ) const
  {
    std::string input_var = section+"/"+var_name;

    if( input.vector_variable_size(input_var) > 1 )
      libmesh_error_msg("ERROR: expression size in input ["+input_var+"] is greater than 1!");
  }

  template<typename Base>
  inline
  void BCFactoryAbstract<Base>::build_var_indices( const MultiphysicsSystem& system,
                                                   const std::vector<std::string>& var_names,
                                                   std::vector<VariableIndex>& var_indices ) const
  {
    var_indices.resize(var_names.size(),libMesh::invalid_uint);

    for( unsigned int i = 0; i < var_names.size(); i++ )
      var_indices[i] = system.variable_number( var_names[i] );
  }

  template<typename Base>
  inline
  void BCFactoryAbstract<Base>::check_state() const
  {
    if( !this->_input )
      libmesh_error_msg("ERROR: must call set_getpot() before building boundary condition!");

    if( !_system )
      libmesh_error_msg("ERROR: must call set_system() before building boundary condition!");

    if( !_bc_ids )
      libmesh_error_msg("ERROR: must call set_bc_ids() before building boundary condition!");

    if( !_fe_var )
      libmesh_error_msg("ERROR: must call set_fe_var() before building boundary condition!");

    if( _section == std::string("DIE!") )
      libmesh_error_msg("ERROR: must call set_section() before building boundary condition!");
  }

  template<typename Base>
  inline
  void BCFactoryAbstract<Base>::reset_state()
  {
    _bc_ids = NULL;
    _fe_var = NULL;
    _section = std::string("DIE!");
  }

  template<typename Base>
  inline
  bool BCFactoryAbstract<Base>::have_bc_type( const std::string& bc_type )
  {
    const std::map<std::string, FactoryAbstract<Base>*>& map =
      FactoryAbstract<Base>::factory_map();

    bool have_bc = false;

    if( map.find(bc_type) != map.end() )
      have_bc = true;

    return have_bc;
  }

} // end namespace GRINS

#endif // GRINS_BC_FACTORY_ABSTRACT_H
