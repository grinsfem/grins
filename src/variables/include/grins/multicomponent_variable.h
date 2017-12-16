//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// GRINS - General Reacting Incompressible Navier-Stokes
//
// Copyright (C) 2014-2017 Paul T. Bauman, Roy H. Stogner
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

#ifndef GRINS_MULTICOMPONENT_VARIABLE_H
#define GRINS_MULTICOMPONENT_VARIABLE_H

// GRINS
#include "grins/fe_variables_base.h"

namespace GRINS
{
  //! Variables that may components, but are not intended to be vectors
  /*! Vector variables should instead use MultcomponentVectorVariable. */
  class MulticomponentVariable : public FEVariablesBase
  {
  public:

    MulticomponentVariable( const std::vector<std::string>& var_names,
                            const std::vector<VariableIndex>& var_indices,
                            const std::set<libMesh::subdomain_id_type>& subdomain_ids )
      : FEVariablesBase(var_names,var_indices,subdomain_ids)
    {}

    virtual ~MulticomponentVariable(){}

    //! Number of components
    unsigned int dim() const
    { return _vars.size(); }

    VariableIndex idx_from_component( unsigned int component ) const
    { libmesh_assert_less( component, _vars.size() );
      return _vars[component]; }

  private:

    MulticomponentVariable();
  };

  class SpeciesMassFractionsVariable : public MulticomponentVariable
  {
  public:

    SpeciesMassFractionsVariable( const std::vector<std::string>& var_names,
                                  const std::vector<VariableIndex>& var_indices,
                                  const std::string& prefix,
                                  const std::string& material,
                                  const std::set<libMesh::subdomain_id_type>& subdomain_ids )
      : MulticomponentVariable(var_names,var_indices,subdomain_ids),
        _prefix(prefix),
        _material(material)
    {}

    virtual ~SpeciesMassFractionsVariable(){}

    unsigned int n_species() const
    { return this->dim(); }

    VariableIndex species( unsigned int species ) const
    { return this->idx_from_component(species); }

    const std::string& material() const
    { return _material; }

    const std::string& prefix() const
    { return _prefix; }

  private:

    SpeciesMassFractionsVariable();

    std::string _prefix;

    std::string _material;

  };

} // end namespace GRINS

#endif // GRINS_MULTICOMPONENT_VARIABLE_H
