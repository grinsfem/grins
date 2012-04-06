//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
// 
// GRINS - General Reacting Incompressible Navier-Stokes 
//
// Copyright (C) 2010-2012 The PECOS Development Team
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the Version 2 GNU General
// Public License as published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301 USA
//
//-----------------------------------------------------------------------el-
//
// $Id$
//
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

#include "solver_factory.h"

GRINS::SolverFactory::SolverFactory(const GetPot& input)
  :_transient( input("unsteady-solver/transient", false ) ),
   _input( input )
{
  return;
}

GRINS::SolverFactory::~SolverFactory()
{
  return;
}

std::tr1::shared_ptr<GRINS::Solver> GRINS::SolverFactory::build()
{
  GRINS::Solver* solver;

  if(_transient)
    {
      solver = new GRINS::UnsteadySolver( _input );
    }
  else
    {
      solver = new GRINS::SteadySolver( _input );
    }

  return std::tr1::shared_ptr<GRINS::Solver>(solver);
}
