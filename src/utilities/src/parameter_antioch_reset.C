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

// Local Includes -----------------------------------
#include "grins/parameter_antioch_reset.h"

#ifdef GRINS_HAVE_ANTIOCH

// Antioch ------------------------------------------
#include "antioch/kinetics_parsing.h"

namespace GRINS
{

ParameterAntiochReset::ParameterAntiochReset() {}


ParameterAntiochReset::ParameterAntiochReset
  (Antioch::ReactionSet<libMesh::Real> & reaction_set,
   const std::string & param_name) :
  _reaction_sets(1, &reaction_set)
{
  std::stringstream stream(param_name);
  std::string keyword;

  std::getline(stream, keyword, '/');
  libmesh_assert(!stream.fail());
  libmesh_assert_equal_to(keyword, "Antioch");

  std::getline(stream, keyword, '/');
  libmesh_assert(!stream.fail());
  _reaction_id = keyword;

  while (std::getline(stream, keyword, '/')) {
    _keywords.push_back(keyword);
  }

  libmesh_assert_greater (_keywords.size(), 0);
}



  /**
   * Setter: change the value of the parameter we access.
   */
void ParameterAntiochReset::set (const libMesh::Number & new_value)
{
  libmesh_assert(!_reaction_sets.empty());
#ifndef NDEBUG
  libMesh::Number val =
    _reaction_sets[0]->get_parameter_of_reaction
      (_reaction_id, _keywords);
#endif
  for (unsigned int i=0; i != _reaction_sets.size(); ++i)
    {
      // If you're already using inconsistent parameters we can't
      // help you.
      libmesh_assert_equal_to
        (val,
         _reaction_sets[i]->get_parameter_of_reaction
           (_reaction_id, _keywords));
      _reaction_sets[i]->set_parameter_of_reaction
        (_reaction_id, _keywords, new_value);
    }
}



/**
 * Getter: get the value of the parameter we access.
 */
const libMesh::Number & ParameterAntiochReset::get () const
  {
    libmesh_assert(!_reaction_sets.empty());
    _current_val =
      _reaction_sets[0]->get_parameter_of_reaction
        (_reaction_id, _keywords);
#ifndef NDEBUG
    // If you're already using inconsistent parameters we can't help
    // you.
    for (unsigned int i=1; i < _reaction_sets.size(); ++i)
      libmesh_assert_equal_to
        (_current_val, _reaction_sets[i]->get_parameter_of_reaction
                         (_reaction_id, _keywords));
#endif
    return _current_val;
  }

} // namespace GRINS

#endif // GRINS_HAVE_ANTIOCH
