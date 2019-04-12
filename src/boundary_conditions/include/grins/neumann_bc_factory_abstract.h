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

#ifndef GRINS_NEUMANN_BC_FACTORY_ABSTRACT_H
#define GRINS_NEUMANN_BC_FACTORY_ABSTRACT_H

// GRINS
#include "grins/bc_factory_abstract.h"
#include "grins/neumann_bc_container.h"

namespace GRINS
{
  // According to the standard, we need declarations of all
  // specializations which precede any automatic instantiations.
  template<> const GetPot* FactoryWithGetPot<NeumannBCContainer>::_input;
  template<> MultiphysicsSystem* BCFactoryAbstract<NeumannBCContainer>::_system;
  template<> const std::set<BoundaryID>* BCFactoryAbstract<NeumannBCContainer>::_bc_ids;
  template<> const FEVariablesBase* BCFactoryAbstract<NeumannBCContainer>::_fe_var;
  template<> std::string BCFactoryAbstract<NeumannBCContainer>::_section;

  class NeumannBCFactoryAbstract : public BCFactoryAbstract<NeumannBCContainer>
  {
  public:
    NeumannBCFactoryAbstract( const std::string& bc_type_name )
      : BCFactoryAbstract<NeumannBCContainer>(bc_type_name),
      _is_homogeneous(false)
    {}

    ~NeumannBCFactoryAbstract(){};

    //! Creates NeumannBCContainer for this Factory object
    /*! This method will handle the creation and population of all
      aspects of the NeumannBCContainer except the function object
      that is used to evaluate the flux, the NeumannBCAbstract
      object. The construction of NeumannBCAbstract subclasses
      is deferred to subclasses of this factory and should
      be implemented in the build_neumann_func method.

      Note that an empty std::unique_ptr<NeumannBCContainer>
      may be returned from create() if the parsed boundary condition
      type is a homogeneous one.  This is allowed since we force
      the user to specify boundary conditions for every
      Variable, for every boundary in the hopes of reducing input
      file errors at runtime. */
    virtual std::unique_ptr<NeumannBCContainer> create();

  protected:

    //! Track if this is a homogeneous Neumann boundary condition
    /*! If so, then we literally do nothing: don't create a NeumannBCAbstract
      object, don't create a NeumannBCContainer, etc. Default is false, so
      subclasses have to opt-in appropriately. */
    bool _is_homogeneous;

    //! Builds the NeumannBCAbstract object for Neumann boundary conditions
    /*! Subclasses should override this function to build the
      NeumannBCAbstract object that corresponds to the variables passed
      in var_names. The FEVariableBase object corresponds to the variable
      associated with this boundary condition, e.g. Velocity. The section
      arguments corresponds to the section to parse for the
      flux in the input file, e.g. input(section+"/"+flux).
      Note that for Variables with more than one component, the flux
      input should be a vector the size of the var_names object, even
      if one of the components is zero. */
    virtual std::shared_ptr<NeumannBCAbstract>
    build_neumann_func( const GetPot& input,
                        MultiphysicsSystem& system,
                        const FEVariablesBase& fe_var,
                        const std::string& section ) =0;

    //! Checks that the flux variable has been set
    /*! Check for both the presence of [section/flux] and that the size
      is the same as var_names. For example, if we are setting a traction vector
      for Displacement variables, then the input flux vector (traction) should
      have the name number of components as var_names. */
    void check_for_flux( const GetPot& input, const std::string& section,
                         const std::vector<std::string>& var_names );

  };

} // end namespace GRINS

#endif // GRINS_NEUMANN_BC_FACTORY_ABSTRACT_H
