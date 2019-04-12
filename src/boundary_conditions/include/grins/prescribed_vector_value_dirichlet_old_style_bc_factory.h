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

#ifndef GRINS_PRESCRIBED_VECTOR_VALUE_DIRICHLET_OLD_STYLE_BC_FACTORY_H
#define GRINS_PRESCRIBED_VECTOR_VALUE_DIRICHLET_OLD_STYLE_BC_FACTORY_H

// GRINS
#include "grins/dirichlet_bc_factory_function_old_style_base.h"

namespace GRINS
{
  // Forward declarations
  class SpeciesMassFractionsVariable;

  class PrescribedVectorValueDirichletOldStyleBCFactory : public DirichletBCFactoryFunctionOldStyleBase<libMesh::FunctionBase<libMesh::Number> >
  {
  public:

    PrescribedVectorValueDirichletOldStyleBCFactory( const std::string& bc_type_name )
      : DirichletBCFactoryFunctionOldStyleBase<libMesh::FunctionBase<libMesh::Number> >(bc_type_name)
    {}

    ~PrescribedVectorValueDirichletOldStyleBCFactory(){};

  protected:

    virtual std::string var_input_string() =0;

    virtual std::unique_ptr<libMesh::FunctionBase<libMesh::Number> >
    build_func( const GetPot& input,
                MultiphysicsSystem& system,
                std::vector<std::string>& var_names,
                const std::string& section );

    virtual void add_funcs( const GetPot& input,
                            MultiphysicsSystem& system,
                            const std::string& input_string,
                            const std::vector<std::string>& var_names,
                            libMesh::CompositeFunction<libMesh::Number>& composite_func ) const;

  };

  class PrescribedVelDirichletOldStyleBCFactory : public PrescribedVectorValueDirichletOldStyleBCFactory
  {
  public:

    PrescribedVelDirichletOldStyleBCFactory( const std::string& bc_type_name )
      : PrescribedVectorValueDirichletOldStyleBCFactory(bc_type_name)
    {}

  protected:

    virtual std::string var_input_string()
    { return "bound_vel"; }

  };

  class PrescribedDispDirichletOldStyleBCFactory : public PrescribedVectorValueDirichletOldStyleBCFactory
  {
  public:

    PrescribedDispDirichletOldStyleBCFactory( const std::string& bc_type_name )
      : PrescribedVectorValueDirichletOldStyleBCFactory(bc_type_name)
    {}

  protected:

    virtual std::string var_input_string()
    { return "displacement"; }

  };

  class PrescribedSpeciesDirichletOldStyleBCFactory : public PrescribedVectorValueDirichletOldStyleBCFactory
  {
  public:

    PrescribedSpeciesDirichletOldStyleBCFactory( const std::string& bc_type_name )
      : PrescribedVectorValueDirichletOldStyleBCFactory(bc_type_name)
    {}

  protected:

    virtual std::string var_input_string()
    { return "bound_species"; }

  };

  class PrescribedMoleFractionsDirichletOldStyleBCFactory : public PrescribedSpeciesDirichletOldStyleBCFactory
  {
  public:

    PrescribedMoleFractionsDirichletOldStyleBCFactory( const std::string& bc_type_name )
      : PrescribedSpeciesDirichletOldStyleBCFactory(bc_type_name)
    {}

  protected:

    virtual void add_funcs( const GetPot& input,
                            MultiphysicsSystem& system,
                            const std::string& input_string,
                            const std::vector<std::string>& var_names,
                            libMesh::CompositeFunction<libMesh::Number>& composite_func ) const;

    template<typename ChemistryType>
    void convert_mole_fracs_and_add_to_func(const GetPot& input,
                                            const std::vector<libMesh::Number>& species_mole_fracs,
                                            const SpeciesMassFractionsVariable& species_fe_var,
                                            libMesh::CompositeFunction<libMesh::Number>& composite_func) const;

  };

} // end namespace GRINS

#endif // GRINS_PRESCRIBED_VECTOR_VALUE_DIRICHLET_OLD_STYLE_BC_FACTORY_H
