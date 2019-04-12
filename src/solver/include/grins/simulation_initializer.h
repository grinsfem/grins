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

#ifndef GRINS_SIMULATION_INITIALIZER_H
#define GRINS_SIMULATION_INITIALIZER_H

namespace GRINS
{
  //! Initialize static objects needed for simulation
  /*! The factory pattern used in GRINS uses static factory objects to handle
    construction of various modules while providing an easy mechanism for
    user extensibility for new variants. However, static linking will strip
    the symbols unless the object is explicitly instantiated. So, we encapluate
    instantiation of all the factory initialization objects here to ensure they
    are instantiated before they are used by Simulation at construction time. */
  class SimulationInitializer
  {
  public:
    SimulationInitializer();
    ~SimulationInitializer(){}

  private:

    // Flag to help guard against multiple initializations
    /*! Simulation objects may need to be created multiple times within a single
      GRINS-linked program (e.g. for parameter variations), so we need to cache
      whether or not this initializer has been called during the current
      program run. If so, then we don't instantiate the factory initializer objects
      because they've already been created. */
    static bool _is_initialized;

  };
} // end namespace GRINS

#endif // GRINS_SIMULATION_INITIALIZER_H
