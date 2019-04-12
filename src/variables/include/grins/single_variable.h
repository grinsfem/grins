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

#ifndef GRINS_SINGLE_VARIABLE_H
#define GRINS_SINGLE_VARIABLE_H

// GRINS
#include "grins/fe_variables_base.h"

namespace GRINS
{
  //! Variables with a single component
  /*! Subclasses can add syntatic sugar, but this class can
    provide the sanity check on there actually being
    one component. */
  class SingleVariable : public FEVariablesBase
  {
  public:

    SingleVariable( const std::vector<std::string>& var_names,
                    const std::vector<VariableIndex>& var_indices,
                    const std::set<libMesh::subdomain_id_type>& subdomain_ids )
      : FEVariablesBase(var_names,var_indices, subdomain_ids )
    {
      if( var_names.size() != 1 )
        libmesh_error_msg("ERROR: SingleVariable must only have a single name!");

      if(var_indices.size() != 1)
        libmesh_error_msg("ERROR: SingleVariable must only have a single index!");
    }

    ~SingleVariable(){}

    VariableIndex var() const
    { libmesh_assert_equal_to( _vars.size(), 1 );
      return this->_vars[0]; }

  private:

    SingleVariable();
  };

  class PrimitiveTempFEVariables : public SingleVariable
  {
  public:
    PrimitiveTempFEVariables( const std::vector<std::string>& var_names,
                              const std::vector<VariableIndex>& var_indices,
                              const std::set<libMesh::subdomain_id_type>& subdomain_ids )
      : SingleVariable(var_names,var_indices,subdomain_ids)
    {}

    ~PrimitiveTempFEVariables(){}

    VariableIndex T() const
    { return this->var(); }

  private:

    PrimitiveTempFEVariables();
  };

  class TurbulenceFEVariables : public SingleVariable
  {
  public:
    TurbulenceFEVariables( const std::vector<std::string>& var_names,
                           const std::vector<VariableIndex>& var_indices,
                           const std::set<libMesh::subdomain_id_type>& subdomain_ids )
      : SingleVariable(var_names,var_indices,subdomain_ids)
    {}

    ~TurbulenceFEVariables(){}

    VariableIndex nu() const
    { return this->var(); }

  private:

    TurbulenceFEVariables();
  };

  class PressureFEVariable : public SingleVariable
  {
  public:
    PressureFEVariable( const std::vector<std::string>& var_names,
                        const std::vector<VariableIndex>& var_indices,
                        const std::set<libMesh::subdomain_id_type>& subdomain_ids )
      : SingleVariable(var_names,var_indices,subdomain_ids)
    {}

    ~PressureFEVariable(){}

    VariableIndex p() const
    { return this->var(); }

  private:

    PressureFEVariable();
  };

  //! Variables with a single SCALAR component
  /*! Will use the type to assert the SCALAR part in the factory. */
  class ScalarVariable : public SingleVariable
  {
  public:
    ScalarVariable( const std::vector<std::string>& var_names,
                    const std::vector<VariableIndex>& var_indices,
                    const std::set<libMesh::subdomain_id_type>& subdomain_ids )
      : SingleVariable(var_names,var_indices,subdomain_ids)
    {}

    ~ScalarVariable(){}

  private:

    ScalarVariable();
  };

  class ThermoPressureVariable : public ScalarVariable
  {
  public:
    ThermoPressureVariable( const std::vector<std::string>& var_names,
                            const std::vector<VariableIndex>& var_indices,
                            const std::set<libMesh::subdomain_id_type>& subdomain_ids )
      : ScalarVariable(var_names,var_indices,subdomain_ids)
    {}

    ~ThermoPressureVariable(){}

    VariableIndex p0() const
    { return this->var(); }

  private:

    ThermoPressureVariable();
  };
} // end namespace GRINS

#endif // GRINS_SINGLE_VARIABLE_H
