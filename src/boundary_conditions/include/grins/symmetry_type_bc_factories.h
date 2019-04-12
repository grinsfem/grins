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

#ifndef GRINS_SYMMETRY_TYPE_BC_FACTORIES_H
#define GRINS_SYMMETRY_TYPE_BC_FACTORIES_H

// GRINS
#include "grins/dirichlet_bc_factory_function_base.h"

// libMesh
#include "libmesh/zero_function.h"

namespace GRINS
{
  class SymmetryTypeBCFactories : public DirichletBCFactoryFunctionBase<libMesh::FunctionBase<libMesh::Number> >
  {
  public:

    SymmetryTypeBCFactories( const std::string& bc_type_name )
      : DirichletBCFactoryFunctionBase(bc_type_name)
    {}

    ~SymmetryTypeBCFactories(){};

  protected:

    //! Trim out the variables that are pinned
    /*! Symmetry-type boundary conditions amount to pinning a subset of
      the group of variables, e.g. Velocity or Displacement. Subclasses
      merely implement this method and remove the variables from var_names
      that *DON'T* need to be pinned.*/
    virtual void trim_var_names( std::vector<std::string>& var_names ) =0;

    //! All the variables are 0, so just return 0 function.
    virtual std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
    build_func( const GetPot& /*input*/,
                MultiphysicsSystem& /*system*/,
                std::vector<std::string>& var_names,
                const std::string& /*section*/ )
    {
      this->trim_var_names(var_names);
      return std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >( new libMesh::ZeroFunction<libMesh::Number> );
    }

  };

  //! Pins x-component of variable (symmetry about yz-plane)
  class YZSymmetryBCFactory : SymmetryTypeBCFactories
  {
  public:
    YZSymmetryBCFactory( const std::string& bc_type_name )
      : SymmetryTypeBCFactories(bc_type_name)
    {}

  protected:
    virtual void trim_var_names( std::vector<std::string>& var_names )
    {
      if( var_names.size() !=3 )
        libmesh_error_msg("ERROR: YZSymmetry requires 3 components in the Variable!");

      // We want to pin the x-component, so remove the last two components
      var_names.pop_back();
      var_names.pop_back();

      libmesh_assert_equal_to(var_names.size(), 1);
    }
  };

  //! Pins y-component of variable (symmetry about xz-plane)
  class XZSymmetryBCFactory : SymmetryTypeBCFactories
  {
  public:
    XZSymmetryBCFactory( const std::string& bc_type_name )
      : SymmetryTypeBCFactories(bc_type_name)
    {}

  protected:
    virtual void trim_var_names( std::vector<std::string>& var_names )
    {
      if( var_names.size() < 2 )
        libmesh_error_msg("ERROR: XZSymmetry requires at least 2 components in the Variable!");

      // We want to pin the y-component, so remove the z-component if it's there
      if( var_names.size() == 3 )
        var_names.pop_back();

      // Now remove x-component
      var_names.erase(var_names.begin());

      libmesh_assert_equal_to(var_names.size(), 1);
    }
  };

  //! Pins z-component of variable (symmetry about xy-plane)
  class XYSymmetryBCFactory : SymmetryTypeBCFactories
  {
  public:
    XYSymmetryBCFactory( const std::string& bc_type_name )
      : SymmetryTypeBCFactories(bc_type_name)
    {}

  protected:
    virtual void trim_var_names( std::vector<std::string>& var_names )
    {
      if( var_names.size() !=3 )
        libmesh_error_msg("ERROR: XYSymmetry requires 3 components in the Variable!");

      // We want to pin the z-component, so remove the x- and y-components
      var_names.erase(var_names.begin());
      var_names.erase(var_names.begin());

      libmesh_assert_equal_to(var_names.size(), 1);
    }
  };

  //! Pins y- and z-component of variable, so can "roll" in the x-direction
  class RollerXBCFactory : SymmetryTypeBCFactories
  {
  public:
    RollerXBCFactory( const std::string& bc_type_name )
      : SymmetryTypeBCFactories(bc_type_name)
    {}

  protected:
    virtual void trim_var_names( std::vector<std::string>& var_names )
    {
      if( var_names.size() < 2 )
        libmesh_error_msg("ERROR: RollerX requires at least 2 components in the Variable!");

      // We want to pin the y- and z-components, so remove the x-component
      var_names.erase(var_names.begin());

      libmesh_assert_less_equal(var_names.size(), 2);
    }
  };

  //! Pins x- and z-component of variable, so can "roll" in the y-direction
  class RollerYBCFactory : SymmetryTypeBCFactories
  {
  public:
    RollerYBCFactory( const std::string& bc_type_name )
      : SymmetryTypeBCFactories(bc_type_name)
    {}

  protected:
    virtual void trim_var_names( std::vector<std::string>& var_names )
    {
      if( var_names.size() < 2 )
        libmesh_error_msg("ERROR: RollerY requires at least 2 components in the Variable!");

      // We want to pin the x- and z-components, so remove the y-component
      var_names.erase(var_names.begin()+1);

      libmesh_assert_less_equal(var_names.size(), 2);
    }
  };

  //! Pins x- and y-component of variable, so can "roll" in the z-direction
  class RollerZBCFactory : SymmetryTypeBCFactories
  {
  public:
    RollerZBCFactory( const std::string& bc_type_name )
      : SymmetryTypeBCFactories(bc_type_name)
    {}

  protected:
    virtual void trim_var_names( std::vector<std::string>& var_names )
    {
      if( var_names.size() !=3 )
        libmesh_error_msg("ERROR: RollerZ requires 3 components in the Variable!");

      // We want to pin the x- and y-components, so remove the z-component
      var_names.pop_back();

      libmesh_assert_equal_to(var_names.size(), 2);
    }
  };

  /*! r-coordinate is assumed to be the x-coordinate in the axisymmetric case,
    so we need to pin the x-coordinate. */
  class AxisymmetryBCFactory : SymmetryTypeBCFactories
  {
  public:
    AxisymmetryBCFactory( const std::string& bc_type_name )
      : SymmetryTypeBCFactories(bc_type_name)
    {}

  protected:
    virtual void trim_var_names( std::vector<std::string>& var_names )
    {
      if( var_names.size() !=2 )
        libmesh_error_msg("ERROR: Axisymmetry requires 2 components in the Variable!");

      // We want to pin the x-components, so remove the y-component
      var_names.pop_back();

      libmesh_assert_equal_to(var_names.size(), 1);
    }
  };

} // end namespace GRINS

#endif // GRINS_SYMMETRY_TYPE_BC_FACTORIES_H
