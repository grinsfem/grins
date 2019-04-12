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

#ifndef GRINS_MULTICOMPONENT_VECTOR_VARIABLE_H
#define GRINS_MULTICOMPONENT_VECTOR_VARIABLE_H

// GRINS
#include "grins/fe_variables_base.h"

namespace GRINS
{
  //! Variables that are effectively vectors
  /*! The variables are effectively vector-valued, but we treat each component
    as a separate variable. This is in contrast to vector-valued FE types,
    e.g. LAGRANGE_VEC, or NEDELEC_ONE.  */
  class MultcomponentVectorVariable : public FEVariablesBase
  {
  public:

    MultcomponentVectorVariable( const std::vector<std::string>& var_names,
                                 const std::vector<VariableIndex>& var_indices,
                                 const std::set<libMesh::subdomain_id_type>& subdomain_ids )
      : FEVariablesBase(var_names,var_indices,subdomain_ids),
        _u_idx(0),
        _v_idx(1),
        _w_idx(2)
    {}

    ~MultcomponentVectorVariable(){}

    //! Number of components
    unsigned int dim() const
    { return _vars.size(); }

    VariableIndex u() const
    { return this->_vars[_u_idx]; }

    VariableIndex v() const
    { return this->_vars[_v_idx]; }

    VariableIndex w() const
    { return this->_vars[_w_idx]; }

  protected:

    unsigned int _u_idx, _v_idx, _w_idx;

  private:

    MultcomponentVectorVariable();
  };

  class DisplacementVariable : public MultcomponentVectorVariable
  {
  public:

    DisplacementVariable( const std::vector<std::string>& var_names,
                          const std::vector<VariableIndex>& var_indices,
                          const std::set<libMesh::subdomain_id_type>& subdomain_ids )
      : MultcomponentVectorVariable(var_names,var_indices,subdomain_ids)
    {}

    ~DisplacementVariable(){}
  };

  class VelocityVariable : public MultcomponentVectorVariable
  {
  public:

    VelocityVariable( const std::vector<std::string>& var_names,
                      const std::vector<VariableIndex>& var_indices,
                      const std::set<libMesh::subdomain_id_type>& subdomain_ids )
      : MultcomponentVectorVariable(var_names,var_indices,subdomain_ids)
    {}

    ~VelocityVariable(){}
  };

} // end namespace GRINS

#endif // GRINS_MULTICOMPONENT_VECTOR_VARIABLE_H
