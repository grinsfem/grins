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

#ifndef GRINS_PARAMETER_ANTIOCH_RESET_H
#define GRINS_PARAMETER_ANTIOCH_RESET_H

#include "grins_config.h"

#ifdef GRINS_HAVE_ANTIOCH

// Antioch ------------------------------------------
#include "antioch/reaction_set.h"

// libMesh ------------------------------------------
#include "libmesh/libmesh_common.h"
#include "libmesh/parameter_accessor.h"

namespace GRINS
{

/**
 * Accessor object allowing reading and modification of Antioch
 * Kinetics variables in a parameter sensitivity calculation.
 */
class ParameterAntiochReset :
  public libMesh::ParameterAccessor<libMesh::Number>
{
public:
  /**
   * Constructor: no parameters attached yet
   */
  ParameterAntiochReset();

  /**
   * Constructor: take the first raw pointer to the parameter
   */
  ParameterAntiochReset
    (Antioch::ReactionSet<libMesh::Real> & reaction_set,
     const std::string & param_name);

  /**
   * A simple reseater won't work with a getter/setter
   */
  virtual ParameterAccessor<libMesh::Number> &
  operator= (libMesh::Number * /* new_ptr */) { libmesh_error(); return *this; }

  /**
   * Setter: change the value of the parameter we access.
   */
  virtual void set (const libMesh::Number & new_value);

  /**
   * Getter: get the value of the parameter we access.
   */
  virtual const libMesh::Number& get () const;

  void push_back
    (Antioch::ReactionSet<libMesh::Real> & new_reaction_set)
    { _reaction_sets.push_back(&new_reaction_set); }

  /**
   * Returns the number of data associated with this parameter.
   * Useful for testing if the resetter is empty/invalid.
   */
  std::size_t size() const { return _reaction_sets.size(); }

  /**
   * Returns a new copy of the accessor.
   */
  virtual libMesh::UniquePtr<ParameterAccessor<libMesh::Number> > clone() const {
    // Default shallow copy works for this class
    ParameterAntiochReset *par = new ParameterAntiochReset(*this);

    return libMesh::UniquePtr<ParameterAccessor<libMesh::Number> >(par);
  }

private:
  // The Antioch getter/setter APIs take a reaction ID and a set of
  // keywords to define which reaction parameter is of concern
  std::string _reaction_id;
  std::vector<std::string> _keywords;

  // We might have multiple Physics, each with its own ReactionSet, in
  // need of perturbation
  std::vector<Antioch::ReactionSet<libMesh::Number> *> _reaction_sets;

  // We need to return a reference from get().  That's a pointless
  // pessimization for libMesh::Number but it might become worthwhile
  // later when we handle field parameters.
  mutable libMesh::Number _current_val;
};

} // namespace GRINS

#endif // GRINS_HAVE_ANTIOCH

#endif // GRINS_PARAMETER_ANTIOCH_RESET_H
